package com.vbuser.apa;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class APAtrapPipeline {

    private static final String CONDA_ENV = "common_env";
    static boolean ignored;

    // 目录与文件
    private static File baseDir;
    private static File scriptsApaDir;
    private static File apaDir;
    private static File bedgraphDir;
    private static File utrDir;
    private static File apaResultsDir;
    private static File refGenomeFile;      // hg38.fa
    private static File chromSizesFile;     // hg38.chrom.sizes
    private static File geneModelBed;        // genemodel.bed

    // APAtrap 工具（Java 实现，不再使用 Perl 可执行文件）
    // 直接调用同包下的 IdentifyDistal3UTR 和 PredictAPA 类

    public static void run(File baseDirectory) throws IOException, InterruptedException {
        baseDir = baseDirectory;
        initDirsAndFiles();

        // 1. （可选）验证 Java 类可用性，已由编译器保证
        System.out.println("使用 Java 重写的 APAtrap 工具。");

        // 2. 生成染色体大小文件
        generateChromSizes();

        // 3. 为所有 BAM 文件生成 bedgraph
        List<File> bedgraphFiles = generateBedgraphsFromBams();

        // 4. 生成过滤后的基因模型 BED 文件（包含单外显子 MSTRG 安全性测试）
        generateFilteredGeneModelBed();

        // 5. 运行 identifyDistal3UTR（Java 版本）
        File refinedUtrBed = runIdentifyDistal3UTR(bedgraphFiles);

        // 6. 运行 predictAPA（Java 高性能实现）
        File apaResultsTxt = runPredictAPA(bedgraphFiles, refinedUtrBed);

        // 7. 运行 Python 下游分析脚本
        runPythonDownstreamAnalysis(apaResultsTxt, refinedUtrBed);
    }

    private static void initDirsAndFiles() throws IOException {
        scriptsApaDir = new File(baseDir, "scripts/apa");
        apaDir = new File(baseDir, "apa");
        bedgraphDir = new File(apaDir, "bedgraph");
        utrDir = new File(apaDir, "utr");
        apaResultsDir = new File(apaDir, "results");

        for (File dir : new File[]{scriptsApaDir, apaDir, bedgraphDir, utrDir, apaResultsDir}) {
            if (!dir.exists() && !dir.mkdirs()) {
                throw new IOException("无法创建目录: " + dir.getAbsolutePath());
            }
        }

        refGenomeFile = new File(baseDir, "hg38.fa");
        if (!refGenomeFile.exists()) {
            throw new IOException("参考基因组不存在: " + refGenomeFile.getAbsolutePath());
        }

        chromSizesFile = new File(baseDir, "hg38.chrom.sizes");
        // gencode.v45.enhanced.gtf
        File enhancedGtf = new File(baseDir, "as/gencode.v45.enhanced.gtf");
        if (!enhancedGtf.exists()) {
            throw new IOException("增强注释文件不存在: " + enhancedGtf.getAbsolutePath());
        }

        geneModelBed = new File(apaDir, "genemodel.bed");

        // 不再需要 Perl 可执行文件路径
    }

    private static void generateChromSizes() throws IOException {
        if (chromSizesFile.exists()) {
            System.out.println("染色体大小文件已存在: " + chromSizesFile);
            return;
        }

        System.out.println("从 FASTA 文件生成染色体大小文件...");
        try (BufferedReader reader = new BufferedReader(new FileReader(refGenomeFile));
             BufferedWriter writer = new BufferedWriter(new FileWriter(chromSizesFile))) {

            String line;
            String currentChrom = null;
            int currentLength = 0;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (currentChrom != null) {
                        writer.write(currentChrom + "\t" + currentLength + "\n");
                    }
                    String header = line.substring(1).trim();
                    int spaceIdx = header.indexOf(' ');
                    currentChrom = (spaceIdx > 0) ? header.substring(0, spaceIdx) : header;
                    currentLength = 0;
                } else {
                    currentLength += line.trim().length();
                }
            }
            if (currentChrom != null) {
                writer.write(currentChrom + "\t" + currentLength + "\n");
            }
        }

        System.out.println("染色体大小文件生成完成: " + chromSizesFile);
    }

    private static List<File> generateBedgraphsFromBams() throws IOException, InterruptedException {
        File bamDir = new File(baseDir, "bam");
        if (!bamDir.exists() || !bamDir.isDirectory()) {
            throw new IOException("BAM 目录不存在: " + bamDir);
        }

        File[] bamFiles = bamDir.listFiles((dir, name) -> name.endsWith(".bam"));
        if (bamFiles == null || bamFiles.length == 0) {
            throw new IOException("在 " + bamDir + " 中未找到任何 BAM 文件");
        }

        List<File> bedgraphFiles = new ArrayList<>();
        for (File bam : bamFiles) {
            String sampleName = bam.getName().replace(".bam", "");
            File bgFile = new File(bedgraphDir, sampleName + ".bedgraph");
            bedgraphFiles.add(bgFile);

            if (bgFile.exists()) {
                System.out.println("BedGraph 已存在，跳过: " + bgFile);
                continue;
            }

            System.out.println("生成 BedGraph: " + sampleName);
            ProcessBuilder pb = new ProcessBuilder(
                    "bedtools", "genomecov",
                    "-bg",
                    "-ibam", bam.getAbsolutePath(),
                    "-g", chromSizesFile.getAbsolutePath(),
                    "-split"
            );
            pb.redirectOutput(bgFile);
            pb.redirectError(ProcessBuilder.Redirect.PIPE);
            Process process = pb.start();

            try (BufferedReader errReader = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
                String line;
                while ((line = errReader.readLine()) != null) {
                    System.err.println("[bedtools] " + line);
                }
            }

            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new IOException("bedtools genomecov 失败，退出码: " + exitCode);
            }
        }
        return bedgraphFiles;
    }

    // ---------- 单外显子 MSTRG 过滤逻辑 ----------
    private static class ExonRecord {
        int start;
        int end;
        ExonRecord(int start, int end) { this.start = start; this.end = end; }
    }

    private static void generateFilteredGeneModelBed() throws IOException, InterruptedException {
        File mergedGtf = new File(baseDir, "as/gencode.v45.enhanced.gtf");
        if (!mergedGtf.exists()) {
            throw new IOException("合并 GTF 不存在: " + mergedGtf);
        }

        Map<String, List<ExonRecord>> transToExons = new LinkedHashMap<>();
        Map<String, String> transToGene = new LinkedHashMap<>();
        Map<String, String> transToChrom = new LinkedHashMap<>();
        Map<String, String> transToStrand = new LinkedHashMap<>();
        parseGtfForTranscripts(mergedGtf, transToExons, transToGene, transToChrom, transToStrand);

        // 分类转录本并过滤单外显子 MSTRG（基于链方向）
        List<String> enstTrans = new ArrayList<>();
        List<String> multiExonMstrg = new ArrayList<>();
        List<String> safeSingleExonMstrg = new ArrayList<>();
        List<String> filteredSingleExonMstrg = new ArrayList<>();

        for (String transId : transToExons.keySet()) {
            List<ExonRecord> exons = transToExons.get(transId);
            String strand = transToStrand.get(transId);
            if (transId.startsWith("ENST")) {
                enstTrans.add(transId);
            } else if (transId.startsWith("MSTRG")) {
                if (exons.size() == 1) {
                    // 单外显子 MSTRG：仅保留明确链方向的转录本
                    if ("+".equals(strand) || "-".equals(strand)) {
                        safeSingleExonMstrg.add(transId);
                    } else {
                        filteredSingleExonMstrg.add(transId);
                    }
                } else {
                    multiExonMstrg.add(transId);
                }
            } else {
                enstTrans.add(transId); // 其他 ID 视为已知转录本
            }
        }

        System.out.println("转录本统计:");
        System.out.println("  ENST (含其他) : " + enstTrans.size());
        System.out.println("  MSTRG 多外显子: " + multiExonMstrg.size());
        System.out.println("  MSTRG 单外显子 (链明确): " + safeSingleExonMstrg.size());
        System.out.println("  MSTRG 单外显子 (链未知，已过滤): " + filteredSingleExonMstrg.size());

        if (!filteredSingleExonMstrg.isEmpty() && safeSingleExonMstrg.isEmpty()) {
            System.out.println("警告: 所有单外显子 MSTRG 链方向未知，将全部舍弃。");
        }

        List<String> finalTransList = new ArrayList<>();
        finalTransList.addAll(enstTrans);
        finalTransList.addAll(multiExonMstrg);
        finalTransList.addAll(safeSingleExonMstrg);

        System.out.println("最终保留转录本总数: " + finalTransList.size());

        File filteredGtf = new File(apaDir, "filtered_enhanced.gtf");
        writeFilteredGtf(mergedGtf, filteredGtf, new HashSet<>(finalTransList));

        generateGeneModelBedFromGtf(filteredGtf);
    }

    private static void parseGtfForTranscripts(File gtfFile,
                                               Map<String, List<ExonRecord>> transToExons,
                                               Map<String, String> transToGene,
                                               Map<String, String> transToChrom,
                                               Map<String, String> transToStrand) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] fields = line.split("\t");
                if (fields.length < 9) continue;
                String feature = fields[2];
                if (!feature.equals("exon")) continue;

                String chrom = fields[0];
                int start = Integer.parseInt(fields[3]);
                int end = Integer.parseInt(fields[4]);
                String strand = fields[6];
                String attr = fields[8];

                String transId = extractAttribute(attr, "transcript_id");
                String geneId = extractAttribute(attr, "gene_id");
                if (transId == null) continue;

                transToExons.computeIfAbsent(transId, k -> new ArrayList<>()).add(new ExonRecord(start, end));
                if (!transToGene.containsKey(transId)) {
                    transToGene.put(transId, geneId != null ? geneId : transId);
                    transToChrom.put(transId, chrom);
                    transToStrand.put(transId, strand);
                }
            }
        }
    }

    private static String extractAttribute(String attr, String key) {
        Pattern p = Pattern.compile(key + " \"([^\"]+)\"");
        Matcher m = p.matcher(attr);
        if (m.find()) {
            return m.group(1);
        }
        return null;
    }

    private static void writeFilteredGtf(File sourceGtf, File destGtf, Set<String> keepTransIds) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(sourceGtf));
             BufferedWriter bw = new BufferedWriter(new FileWriter(destGtf))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) {
                    bw.write(line + "\n");
                    continue;
                }
                String[] fields = line.split("\t");
                if (fields.length < 9) {
                    bw.write(line + "\n");
                    continue;
                }
                String attr = fields[8];
                String transId = extractAttribute(attr, "transcript_id");
                if (transId != null && keepTransIds.contains(transId)) {
                    bw.write(line + "\n");
                } else if (!attr.contains("transcript_id")) {
                    bw.write(line + "\n");
                }
            }
        }
    }

    private static void generateGeneModelBedFromGtf(File gtfFile) throws IOException, InterruptedException {
        generateGeneModelBedFromGtf(gtfFile, geneModelBed);
    }

    private static void generateGeneModelBedFromGtf(File gtfFile, File outputBed) throws IOException, InterruptedException {
        File gtfToGenePred = ensureUcscTool("gtfToGenePred");
        File genePredToBed = ensureUcscTool("genePredToBed");
        File genePredFile = new File(outputBed.getParent(), outputBed.getName().replace(".bed", ".genePred"));

        runCommand(gtfToGenePred.getAbsolutePath(), gtfFile.getAbsolutePath(), genePredFile.getAbsolutePath());
        runCommand(genePredToBed.getAbsolutePath(), genePredFile.getAbsolutePath(), outputBed.getAbsolutePath());
        ignored = genePredFile.delete();
    }

    private static File ensureUcscTool(String toolName) throws IOException, InterruptedException {
        File toolFile = new File(scriptsApaDir, toolName);
        if (toolFile.exists()) {
            ignored = toolFile.setExecutable(true);
            return toolFile;
        }

        String url = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/" + toolName;
        System.out.println("下载 UCSC 工具: " + toolName);
        runCommand("wget", "-O", toolFile.getAbsolutePath(), url);
        ignored = toolFile.setExecutable(true);
        return toolFile;
    }

    // ---------- 调用 Java 版 identifyDistal3UTR ----------
    private static File runIdentifyDistal3UTR(List<File> bedgraphFiles) throws IOException {
        File outputBed = new File(utrDir, "refined_3utr.bed");
        if (outputBed.exists()) {
            System.out.println("refined_3utr.bed 已存在，跳过 identifyDistal3UTR。");
            return outputBed;
        }

        List<String> argsList = new ArrayList<>();
        argsList.add("-i");
        for (File f : bedgraphFiles) {
            argsList.add(f.getAbsolutePath());
        }
        argsList.add("-m");
        argsList.add(geneModelBed.getAbsolutePath());
        argsList.add("-o");
        argsList.add(outputBed.getAbsolutePath());
        // 可根据需要添加其他参数，如 -w, -e, -c, -p, -s

        System.out.println("运行 identifyDistal3UTR (Java 实现)...");
        String[] args = argsList.toArray(new String[0]);
        try {
            IdentifyDistal3UTR.main(args);
        } catch (Exception e) {
            throw new IOException("IdentifyDistal3UTR 执行失败: " + e.getMessage(), e);
        }

        if (!outputBed.exists()) {
            throw new IOException("IdentifyDistal3UTR 未能生成输出文件: " + outputBed);
        }
        return outputBed;
    }

    // ---------- 调用 Java 版 predictAPA ----------
    private static File runPredictAPA(List<File> bedgraphFiles, File utrBed) throws IOException {
        File outputTxt = new File(apaResultsDir, "all_apa_results.txt");
        if (outputTxt.exists()) {
            System.out.println("all_apa_results.txt 已存在，跳过 predictAPA。");
            return outputTxt;
        }

        int numSamples = bedgraphFiles.size();
        List<String> argsList = new ArrayList<>();
        argsList.add("-i");
        for (File f : bedgraphFiles) {
            argsList.add(f.getAbsolutePath());
        }
        argsList.add("-g");
        argsList.add("1");
        argsList.add("-n");
        argsList.add(String.valueOf(numSamples));
        argsList.add("-u");
        argsList.add(utrBed.getAbsolutePath());
        argsList.add("-o");
        argsList.add(outputTxt.getAbsolutePath());

        System.out.println("运行 PredictAPA (Java 高性能实现) ...");
        String[] args = argsList.toArray(new String[0]);
        try {
            PredictAPA.main(args);
        } catch (Exception e) {
            throw new IOException("PredictAPA 执行失败: " + e.getMessage(), e);
        }

        if (!outputTxt.exists()) {
            throw new IOException("PredictAPA 未能生成输出文件: " + outputTxt);
        }
        return outputTxt;
    }

    private static void runPythonDownstreamAnalysis(File apaResults, File utrBed) throws IOException, InterruptedException {
        File predictionsFile = new File(apaResultsDir, "apa_predictions.txt");
        Files.copy(apaResults.toPath(), predictionsFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

        File pythonScript = new File(new File(baseDir, "scripts/apa"), "run_apa_pipeline.py");
        if (!pythonScript.exists()) {
            System.err.println("警告: Python 脚本 run_apa_pipeline.py 不存在，跳过下游分析。");
            return;
        }

        File outDir = new File(new File(baseDir, "apa"),"results");
        if (!outDir.exists() && !outDir.mkdirs()) {
            throw new IOException("无法创建输出目录: " + outDir);
        }

        List<String> cmd = new ArrayList<>();
        cmd.add("conda");
        cmd.add("run");
        cmd.add("--no-capture-output");
        cmd.add("-n");
        cmd.add(CONDA_ENV);
        cmd.add("python");
        cmd.add(pythonScript.getAbsolutePath());
        cmd.add("--all_results"); cmd.add(apaResults.getAbsolutePath());
        cmd.add("--predictions"); cmd.add(predictionsFile.getAbsolutePath());
        cmd.add("--utr"); cmd.add(utrBed.getAbsolutePath());
        cmd.add("--genome"); cmd.add(refGenomeFile.getAbsolutePath());
        cmd.add("--outdir"); cmd.add(outDir.getAbsolutePath());

        System.out.println("执行下游 Python 分析...");
        System.out.println("命令: " + String.join(" ", cmd));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.environment().put("PYTHONUNBUFFERED", "1");
        pb.redirectErrorStream(true);
        Process process = pb.start();

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }
        }

        int exitCode = process.waitFor();
        if (exitCode != 0) {
            throw new IOException("Python 分析脚本执行失败，退出码: " + exitCode);
        }
        System.out.println("APA 下游分析完成，结果保存于: " + outDir);
    }

    // ---------- 通用命令执行 ----------
    private static void runCommand(String... cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        Process process = pb.start();
        String output = readStream(process.getInputStream());
        int exitCode = process.waitFor();
        if (exitCode != 0) {
            throw new IOException("命令执行失败: " + String.join(" ", cmd) + "\n" + output);
        }
    }

    private static String readStream(InputStream is) throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                sb.append(line).append("\n");
            }
            return sb.toString();
        }
    }
}