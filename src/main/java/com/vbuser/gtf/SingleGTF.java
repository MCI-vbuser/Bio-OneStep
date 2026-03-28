package com.vbuser.gtf;

import java.io.*;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.stream.Collectors;

public class SingleGTF {

    // ======================== 配置参数 ========================
    private static final String CONDA_ENV = "common_env";
    private static final String GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz";
    private static final String GENCODE_GTF = "gencode.v45.annotation.gtf";
    private static final String REF_GENOME = "hg38.fa";
    private static final int THREADS = 4;

    // 输出目录
    private static File baseDir;
    private static File sraConvertDir;
    private static File cleanedFastqDir;
    private static File bamDir;
    private static File gtfSingleDir;
    private static File mergeDir;
    private static File quantDir;
    private static File gffcompareDir;
    static boolean ignored;     //这个布尔值是为了避免IDEA对文件操作报Warning

    // 样本信息
    private static final List<SampleInfo> sampleInfos = new ArrayList<>();
    // 存储清洗后的 FastQ 文件，key=样本名，value=文件列表（单端或双端）
    private static final Map<String, List<File>> cleanedFastqMap = new HashMap<>();

    private static class SampleInfo {
        String sampleName;
        File sortedBam;
        File rawGtf;
    }

    // ======================== 初始化 ========================
    private static void initDirs() throws IOException {
        if (baseDir == null) {
            String jarPath = SingleGTF.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            jarPath = URLDecoder.decode(jarPath, StandardCharsets.UTF_8.name());
            File jarFile = new File(jarPath);
            baseDir = jarFile.getParentFile();
            if (baseDir == null) baseDir = new File(".");
        }

        sraConvertDir = new File(baseDir, "sra_fastq");
        cleanedFastqDir = new File(baseDir, "cleaned_fastq");
        bamDir = new File(baseDir, "bam");
        gtfSingleDir = new File(baseDir, "gtf_single");
        mergeDir = new File(baseDir, "gtf_merge");
        quantDir = new File(mergeDir, "quant");
        gffcompareDir = new File(mergeDir, "gffcompare");

        for (File dir : new File[]{sraConvertDir, cleanedFastqDir, bamDir, gtfSingleDir, mergeDir, quantDir, gffcompareDir}) {
            if (!dir.exists() && !dir.mkdirs()) {
                throw new IOException("无法创建目录: " + dir.getAbsolutePath());
            }
        }
    }

    // ======================== 下载 GENCODE ========================
    private static void downloadGencode() throws IOException, InterruptedException {
        File gencode = new File(baseDir, GENCODE_GTF);
        if (!gencode.exists()) {
            System.out.println("下载 GENCODE v45 注释...");
            runCommand("wget", "-O", GENCODE_GTF + ".gz", GENCODE_URL);
            runCommand("gunzip", GENCODE_GTF + ".gz");
            ignored = new File(GENCODE_GTF + ".gz").delete();
            System.out.println("下载完成: " + gencode.getAbsolutePath());
        } else {
            System.out.println("注释文件已存在: " + gencode.getAbsolutePath());
        }
    }

    // ======================== SRA → FastQ ========================
    private static void convertSraToFastq(File sraFile) throws IOException, InterruptedException {
        String baseName = sraFile.getName().replaceAll("\\.sra$", "");
        File out1 = new File(sraConvertDir, baseName + "_1.fastq");
        File out2 = new File(sraConvertDir, baseName + "_2.fastq");
        if (out1.exists() && out2.exists()) {
            System.out.println("FastQ 已存在，跳过转换: " + baseName);
            return;
        }
        System.out.println("转换 SRA: " + sraFile.getName());
        runCommand("fastq-dump", "--split-files", sraFile.getAbsolutePath(), "-O", sraConvertDir.getAbsolutePath());
    }

    // ======================== fastp 清洗 ========================
    private static void cleanFastq(File fq1, File fq2, String sampleName) throws IOException, InterruptedException {
        File out1 = new File(cleanedFastqDir, sampleName + "_1_clean.fastq");
        File out2 = fq2 != null ? new File(cleanedFastqDir, sampleName + "_2_clean.fastq") : null;
        File html = new File(cleanedFastqDir, sampleName + "_fastp_report.html");
        File json = new File(cleanedFastqDir, sampleName + "_fastp_report.json");

        if (out1.exists() && (out2 == null || out2.exists())) {
            System.out.println("清洗后 FastQ 已存在，跳过: " + sampleName);
            // 仍然记录到 map，供后续使用
            cleanedFastqMap.computeIfAbsent(sampleName, k -> new ArrayList<>()).add(out1);
            if (out2 != null) cleanedFastqMap.get(sampleName).add(out2);
            return;
        }

        List<String> cmd = new ArrayList<>();
        cmd.add("fastp");
        cmd.add("-i"); cmd.add(fq1.getAbsolutePath());
        if (fq2 != null) {
            cmd.add("-I"); cmd.add(fq2.getAbsolutePath());
            cmd.add("-o"); cmd.add(out1.getAbsolutePath());
            cmd.add("-O"); cmd.add(out2.getAbsolutePath());
        } else {
            cmd.add("-o"); cmd.add(out1.getAbsolutePath());
        }
        cmd.add("-q"); cmd.add("20");
        cmd.add("-u"); cmd.add("30");
        cmd.add("-l"); cmd.add("50");
        cmd.add("--n_base_limit"); cmd.add("5");
        cmd.add("--trim_poly_g");
        cmd.add("--trim_poly_x");
        cmd.add("--trim_front1"); cmd.add("1");
        cmd.add("--trim_tail1"); cmd.add("1");
        cmd.add("-w"); cmd.add(String.valueOf(THREADS));
        cmd.add("-h"); cmd.add(html.getAbsolutePath());
        cmd.add("-j"); cmd.add(json.getAbsolutePath());

        System.out.println("运行 fastp 清洗: " + sampleName);
        runCommandWithConda(cmd.toArray(new String[0]));

        // 记录清洗后的文件
        cleanedFastqMap.computeIfAbsent(sampleName, k -> new ArrayList<>()).add(out1);
        if (out2 != null) cleanedFastqMap.get(sampleName).add(out2);
    }

    // ======================== minimap2 比对 → sorted BAM ========================
    private static void alignToBam(File fq1, File fq2, String sampleName) throws IOException, InterruptedException {
        File finalBam = new File(bamDir, sampleName + ".bam");
        if (finalBam.exists()) {
            System.out.println("BAM 已存在，跳过比对: " + sampleName);
            SampleInfo info = new SampleInfo();
            info.sampleName = sampleName;
            info.sortedBam = finalBam;
            sampleInfos.add(info);
            return;
        }

        File tmpSam = new File(bamDir, sampleName + ".sam");
        List<String> cmd = new ArrayList<>();
        cmd.add("minimap2");
        cmd.add("-ax"); cmd.add("splice");
        cmd.add("-t"); cmd.add(String.valueOf(THREADS));
        cmd.add(REF_GENOME);
        cmd.add(fq1.getAbsolutePath());
        if (fq2 != null) cmd.add(fq2.getAbsolutePath());
        cmd.add("-o"); cmd.add(tmpSam.getAbsolutePath());

        System.out.println("minimap2 比对: " + sampleName);
        runCommandWithConda(cmd.toArray(new String[0]));

        System.out.println("排序并转换 BAM: " + sampleName);
        runCommandWithConda("samtools", "sort", "-@", String.valueOf(THREADS),
                "-o", finalBam.getAbsolutePath(), tmpSam.getAbsolutePath());
        runCommandWithConda("samtools", "index", finalBam.getAbsolutePath());
        ignored = tmpSam.delete();

        SampleInfo info = new SampleInfo();
        info.sampleName = sampleName;
        info.sortedBam = finalBam;
        sampleInfos.add(info);
    }

    // ======================== 参考引导组装 ========================
    private static void assembleGuided(SampleInfo info) throws IOException, InterruptedException {
        File rawGtf = new File(gtfSingleDir, info.sampleName + ".gtf");
        if (rawGtf.exists()) {
            System.out.println("GTF 已存在，跳过组装: " + info.sampleName);
            info.rawGtf = rawGtf;
            return;
        }

        File refGtf = new File(baseDir, GENCODE_GTF);
        if (!refGtf.exists()) throw new IOException("参考注释缺失: " + refGtf);

        System.out.println("StringTie 参考引导组装: " + info.sampleName);
        runCommandWithConda("stringtie",
                "-o", rawGtf.getAbsolutePath(),
                "-G", refGtf.getAbsolutePath(),
                "-p", String.valueOf(THREADS),
                "-l", "STRG",
                info.sortedBam.getAbsolutePath());

        info.rawGtf = rawGtf;
    }

    // ======================== 合并所有样本的 GTF ========================
    private static void mergeGtfs() throws IOException, InterruptedException {
        if (sampleInfos.isEmpty()) return;
        File mergedGtf = new File(mergeDir, "merged.gtf");
        if (mergedGtf.exists()) {
            System.out.println("合并文件已存在: " + mergedGtf.getAbsolutePath());
            return;
        }

        List<String> cmd = new ArrayList<>();
        cmd.add("stringtie");
        cmd.add("--merge");
        cmd.add("-o");
        cmd.add(mergedGtf.getAbsolutePath());
        for (SampleInfo info : sampleInfos) cmd.add(info.rawGtf.getAbsolutePath());

        System.out.println("合并 GTF...");
        runCommandWithConda(cmd.toArray(new String[0]));
        System.out.println("合并后大小: " + mergedGtf.length() + " 字节");
    }

    // ======================== 量化并运行 gffcompare ========================
    private static void quantifyAndCompare() throws IOException, InterruptedException {
        File mergedGtf = new File(mergeDir, "merged.gtf");
        if (!mergedGtf.exists()) throw new IOException("合并 GTF 不存在");

        for (SampleInfo info : sampleInfos) {
            File sampleQuantDir = new File(quantDir, info.sampleName);
            if (!sampleQuantDir.exists() && !sampleQuantDir.mkdirs())
                throw new IOException("无法创建量化目录: " + sampleQuantDir);

            File ballgownGtf = new File(sampleQuantDir, "ballgown_" + info.sampleName + ".gtf");
            File abundFile = new File(sampleQuantDir, "abund_" + info.sampleName + ".tsv");

            if (!ballgownGtf.exists() || !abundFile.exists()) {
                System.out.println("量化样本: " + info.sampleName);
                runCommandInDir(sampleQuantDir,
                        "conda", "run", "-n", CONDA_ENV,
                        "stringtie", "-e", "-B", "-p", String.valueOf(THREADS),
                        "-G", mergedGtf.getAbsolutePath(),
                        "-o", ballgownGtf.getAbsolutePath(),
                        "-A", abundFile.getAbsolutePath(),
                        info.sortedBam.getAbsolutePath());
            }

            runGffcompare(ballgownGtf, info.sampleName);
        }
    }

    private static void runGffcompare(File ballgownGtf, String sampleName) throws IOException, InterruptedException {
        if (!ballgownGtf.exists() || ballgownGtf.length() == 0) {
            System.err.println("ballgown GTF 为空，跳过 gffcompare: " + sampleName);
            return;
        }

        File sampleGffDir = new File(gffcompareDir, sampleName);
        if (!sampleGffDir.exists() && !sampleGffDir.mkdirs())
            throw new IOException("无法创建 gffcompare 目录");

        File localGtf = new File(sampleGffDir, ballgownGtf.getName());
        if (!localGtf.exists())
            Files.copy(ballgownGtf.toPath(), localGtf.toPath(), StandardCopyOption.REPLACE_EXISTING);

        File refGtf = new File(baseDir, GENCODE_GTF);
        if (!refGtf.exists()) throw new IOException("参考注释缺失");

        String prefix = new File(sampleGffDir, "g_" + sampleName).getAbsolutePath();
        List<String> cmd = Arrays.asList("conda", "run", "-n", CONDA_ENV,
                "gffcompare", "-r", refGtf.getAbsolutePath(),
                "-o", prefix, localGtf.getAbsolutePath());

        System.out.println("执行 gffcompare: " + String.join(" ", cmd));
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.directory(sampleGffDir);
        pb.redirectErrorStream(true);
        Process process = pb.start();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            reader.lines().forEach(System.out::println);
        }
        int exitCode = process.waitFor();
        if (exitCode != 0) throw new IOException("gffcompare 失败，退出码: " + exitCode);
        System.out.println("gffcompare 完成，结果在: " + sampleGffDir);
    }

    // ======================== 主入口 ========================
    public static void handle(File input, boolean sortedEnsured) {
        try {
            initDirs();
            downloadGencode();

            // 收集所有输入文件
            List<File> allFiles = new ArrayList<>();
            if (input.isFile()) allFiles.add(input);
            else if (input.isDirectory()) collectFiles(input, allFiles);
            else throw new IOException("无效输入");

            // 分类
            List<File> sraList = allFiles.stream().filter(f -> f.getName().endsWith(".sra")).collect(Collectors.toList());
            List<File> fastqList = allFiles.stream().filter(f -> f.getName().endsWith(".fastq") || f.getName().endsWith(".fq")).collect(Collectors.toList());
            List<File> bamList = allFiles.stream().filter(f -> f.getName().endsWith(".bam")).collect(Collectors.toList());

            // 1. SRA → FastQ
            for (File sra : sraList) convertSraToFastq(sra);

            // 2. 收集所有原始 FastQ（包括转换后和已有的）
            List<File> allFastq = new ArrayList<>();
            collectFiles(sraConvertDir, allFastq);
            allFastq.addAll(fastqList);
            Map<String, List<File>> rawPairs = groupFastqPairs(allFastq);

            // 3. fastp 清洗（直接填充 cleanedFastqMap）
            for (Map.Entry<String, List<File>> e : rawPairs.entrySet()) {
                List<File> pair = e.getValue();
                File fq1 = pair.get(0);
                File fq2 = pair.size() > 1 ? pair.get(1) : null;
                cleanFastq(fq1, fq2, e.getKey());
            }

            // 4. 使用 cleanedFastqMap 进行比对
            for (Map.Entry<String, List<File>> e : cleanedFastqMap.entrySet()) {
                String sampleName = e.getKey();
                List<File> pair = e.getValue();
                File fq1 = pair.get(0);
                File fq2 = pair.size() > 1 ? pair.get(1) : null;
                alignToBam(fq1, fq2, sampleName);
            }

            // 5. 处理直接输入的 BAM（若有）
            for (File bam : bamList) {
                String sampleName = bam.getName().replaceAll("\\.bam$", "");
                SampleInfo info = new SampleInfo();
                info.sampleName = sampleName;
                if (sortedEnsured) {
                    info.sortedBam = bam;
                } else {
                    File sortedBam = new File(bam.getParent(), sampleName + ".sorted.bam");
                    if (!sortedBam.exists()) {
                        runCommandWithConda("samtools", "sort", "-@", String.valueOf(THREADS),
                                "-o", sortedBam.getAbsolutePath(), bam.getAbsolutePath());
                    }
                    info.sortedBam = sortedBam;
                }
                sampleInfos.add(info);
            }

            // 6. 参考引导组装
            for (SampleInfo info : sampleInfos) assembleGuided(info);

            // 7. 合并
            mergeGtfs();

            // 8. 量化 + gffcompare
            quantifyAndCompare();

        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    // ======================== 辅助方法 ========================
    private static void collectFiles(File dir, List<File> list) {
        if (dir == null) return;
        File[] children = dir.listFiles();
        if (children != null) {
            for (File f : children) {
                if (f.isFile()) list.add(f);
                else if (f.isDirectory()) collectFiles(f, list);
            }
        }
    }

    private static Map<String, List<File>> groupFastqPairs(List<File> fastqFiles) {
        Map<String, List<File>> map = new HashMap<>();
        for (File f : fastqFiles) {
            String name = f.getName();
            // 匹配标准格式：sample_1.fastq 或 sample_2.fastq
            String base = name.replaceAll("(_[12])?\\.(fastq|fq)$", "");
            map.computeIfAbsent(base, k -> new ArrayList<>()).add(f);
        }
        for (List<File> list : map.values()) list.sort(Comparator.comparing(File::getName));
        return map;
    }

    private static void runCommand(String... cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        Process p = pb.start();
        String out = readStream(p.getInputStream());
        int code = p.waitFor();
        if (code != 0) throw new IOException("命令失败: " + String.join(" ", cmd) + "\n" + out);
    }

    private static void runCommandWithConda(String... cmd) throws IOException, InterruptedException {
        List<String> full = new ArrayList<>();
        full.add("conda"); full.add("run"); full.add("-n"); full.add(CONDA_ENV);
        full.addAll(Arrays.asList(cmd));
        runCommand(full.toArray(new String[0]));
    }

    private static void runCommandInDir(File dir, String... cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.directory(dir);
        pb.redirectErrorStream(true);
        Process p = pb.start();
        String out = readStream(p.getInputStream());
        int code = p.waitFor();
        if (code != 0) throw new IOException("命令失败: " + String.join(" ", cmd) + "\n" + out);
    }

    private static String readStream(InputStream is) throws IOException {
        try (BufferedReader r = new BufferedReader(new InputStreamReader(is))) {
            return r.lines().collect(Collectors.joining("\n"));
        }
    }
}