package com.vbuser.gtf;

import java.io.*;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SingleGTF {

    // 与 FastqProcessor 保持一致
    private static final String CONDA_ENV = "common_env";
    private static final String GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz";
    private static final String GENCODE_GTF = "gencode.v45.annotation.gtf";
    static boolean ignored;

    // GTF 输出目录（相对于 JAR 所在目录）
    private static File gtfOutputDir;

    /**
     * 获取 GTF 输出目录（相对于 JAR 文件所在目录的 ./gtf_single）
     * 如果目录不存在则创建
     */
    private static synchronized File getGtfOutputDir() throws IOException {
        if (gtfOutputDir == null) {
            // 获取 JAR 文件所在目录
            String jarPath = SingleGTF.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            // 处理 URL 编码（例如空格 %20）
            jarPath = URLDecoder.decode(jarPath, StandardCharsets.UTF_8.name());
            File jarFile = new File(jarPath);
            File jarDir = jarFile.getParentFile();
            if (jarDir == null) {
                // 如果无法获取，回退到当前工作目录
                jarDir = new File(".");
            }
            gtfOutputDir = new File(jarDir, "gtf_single");
            if (!gtfOutputDir.exists() && !gtfOutputDir.mkdirs()) {
                throw new IOException("无法创建 GTF 输出目录: " + gtfOutputDir.getAbsolutePath());
            }
            System.out.println("GTF 输出目录: " + gtfOutputDir.getAbsolutePath());
        }
        return gtfOutputDir;
    }

    /**
     * 获取合并输出目录（./gtf_merge），若不存在则创建
     */
    private static synchronized File getMergeDir() throws IOException {
        String jarPath = SingleGTF.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        jarPath = URLDecoder.decode(jarPath, StandardCharsets.UTF_8.name());
        File jarFile = new File(jarPath);
        File jarDir = jarFile.getParentFile();
        if (jarDir == null) jarDir = new File(".");
        File mergeDir = new File(jarDir, "gtf_merge");
        if (!mergeDir.exists() && !mergeDir.mkdirs()) {
            throw new IOException("无法创建合并输出目录: " + mergeDir.getAbsolutePath());
        }
        return mergeDir;
    }

    /**
     * 递归收集目录下所有 .bam 文件
     */
    private static void collectBamFiles(File dir, List<File> list) {
        File[] children = dir.listFiles();
        if (children != null) {
            for (File f : children) {
                if (f.isFile() && f.getName().endsWith(".bam")) {
                    list.add(f);
                } else if (f.isDirectory()) {
                    collectBamFiles(f, list);
                }
            }
        }
    }

    /**
     * 下载 Gencode 注释文件（如果不存在）
     */
    private static void download_gencode_annotation() throws IOException, InterruptedException {
        File gencode = new File(GENCODE_GTF);
        if (!gencode.exists()) {
            System.out.println("正在下载 " + GENCODE_GTF);
            runCommand("wget", "-O", GENCODE_GTF + ".gz", GENCODE_URL);
            runCommand("gunzip", GENCODE_GTF + ".gz");
            System.out.println("下载完成: " + GENCODE_GTF);
        } else {
            System.out.println("注释文件已存在: " + GENCODE_GTF);
        }
    }

    /**
     * 处理单个 BAM 文件或递归处理目录
     * 当输入为目录时，先收集所有 BAM，逐个处理，最后若 BAM 数量 ≥ 2 则执行合并
     *
     * @param bam_files      BAM 文件或目录
     * @param sorted_ensured 是否确保 BAM 已排序
     */
    public static void handle(File bam_files, boolean sorted_ensured) {
        System.out.println(sorted_ensured ? "确保 BAM 已排序" : "不确保 BAM 已排序");
        try {
            download_gencode_annotation();
            if (bam_files.isFile() && bam_files.getName().endsWith(".bam")) {
                processSingleBam(bam_files, sorted_ensured);
            } else if (bam_files.isDirectory()) {
                List<File> bamList = new ArrayList<>();
                collectBamFiles(bam_files, bamList);
                for (File bam : bamList) {
                    processSingleBam(bam, sorted_ensured);
                }
                if (bamList.size() >= 2) {
                    mergeAndCompare();
                } else if (bamList.size() == 1) {
                    System.out.println("只有一个 BAM 文件，跳过合并。");
                } else {
                    System.out.println("未找到任何 BAM 文件。");
                }
            }
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * 处理单个 BAM 文件：用 stringtie 组装 GTF，并移动到统一输出目录
     *
     * @param bam            BAM 文件
     * @param sorted_ensured 是否确保排序
     */
    private static void processSingleBam(File bam, boolean sorted_ensured) throws IOException, InterruptedException {
        String bamPath = bam.getAbsolutePath();
        // 目标 GTF 文件名：与 BAM 同名，扩展名改为 .gtf
        String gtfFileName = bam.getName().replaceAll("\\.bam$", ".gtf");
        File targetGtf = new File(getGtfOutputDir(), gtfFileName);

        // 如果 GTF 已存在，跳过处理
        if (targetGtf.exists()) {
            System.out.println("GTF 文件已存在，跳过: " + targetGtf.getName());
            return;
        }

        // 如果需要排序且未确保排序，则先对 BAM 排序
        String inputBam = bamPath;
        if (!sorted_ensured) {
            System.out.println("对 BAM 进行排序: " + bam.getName());
            String sortedBam = bamPath.replaceAll("\\.bam$", ".sorted.bam");
            runCommandWithConda("samtools", "sort", bamPath, "-o", sortedBam);
            inputBam = sortedBam;
        }

        // 运行 stringtie 组装 GTF，直接输出到目标目录
        System.out.println("组装 GTF: " + bam.getName() + " -> " + targetGtf.getAbsolutePath());
        runCommandWithConda("stringtie", "-o", targetGtf.getAbsolutePath(), "-G", GENCODE_GTF, inputBam);

        // 如果生成了临时排序文件，则删除
        if (!sorted_ensured && !inputBam.equals(bamPath)) {
            ignored = new File(inputBam).delete();
            System.out.println("已删除临时排序文件: " + inputBam);
        }

        System.out.println("完成: " + targetGtf.getAbsolutePath());
    }

    /**
     * 合并所有样本的 GTF，并为每个样本生成与合并 GTF 的对应关系
     */
    private static void mergeAndCompare() throws IOException, InterruptedException {
        // 获取 GTF 输出目录（所有样本 GTF 所在位置）
        File gtfDir = getGtfOutputDir();
        File[] gtfFiles = gtfDir.listFiles((dir, name) -> name.endsWith(".gtf"));
        if (gtfFiles == null || gtfFiles.length < 2) {
            System.out.println("GTF 文件数量不足 2 个，跳过合并");
            return;
        }

        // 创建合并输出目录
        File mergeDir = getMergeDir();
        File mergedGtf = new File(mergeDir, "merged.gtf");

        // 1. 使用 stringtie --merge 合并所有 GTF
        System.out.println("正在合并 " + gtfFiles.length + " 个 GTF 文件...");
        List<String> mergeCmd = new ArrayList<>();
        mergeCmd.add("stringtie");
        mergeCmd.add("--merge");
        mergeCmd.add("-o");
        mergeCmd.add(mergedGtf.getAbsolutePath());
        for (File gtf : gtfFiles) {
            mergeCmd.add(gtf.getAbsolutePath());
        }
        runCommandWithConda(mergeCmd.toArray(new String[0]));

        // 2. 对每个样本 GTF 运行 gffcompare
        for (File sampleGtf : gtfFiles) {
            String sampleName = sampleGtf.getName().replaceAll("\\.gtf$", "");
            File sampleOutDir = new File(mergeDir, sampleName);
            if (!sampleOutDir.exists() && !sampleOutDir.mkdirs()) {
                throw new IOException("无法创建样本输出目录: " + sampleOutDir.getAbsolutePath());
            }
            String outputPrefix = new File(sampleOutDir, sampleName).getAbsolutePath();
            System.out.println("正在比较 " + sampleGtf.getName() + " 与合并 GTF...");
            runCommandWithConda("gffcompare", "-r", mergedGtf.getAbsolutePath(),
                    "-o", outputPrefix, sampleGtf.getAbsolutePath());
        }
        System.out.println("合并与比较完成，结果保存在: " + mergeDir.getAbsolutePath());
    }

    /**
     * 执行系统命令（不使用 conda 环境）
     */
    private static void runCommand(String... cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        Process process = pb.start();
        int exitCode = process.waitFor();
        if (exitCode != 0) {
            String error = readStream(process.getInputStream());
            throw new IOException("命令执行失败: " + String.join(" ", cmd) + "\n" + error);
        }
    }

    /**
     * 在 conda 环境中执行命令
     */
    private static void runCommandWithConda(String... cmd) throws IOException, InterruptedException {
        List<String> fullCmd = new ArrayList<>();
        fullCmd.add("conda");
        fullCmd.add("run");
        fullCmd.add("-n");
        fullCmd.add(CONDA_ENV);
        fullCmd.addAll(Arrays.asList(cmd));
        runCommand(fullCmd.toArray(new String[0]));
    }

    /**
     * 读取输入流内容（用于错误输出）
     */
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