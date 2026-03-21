package com.vbuser.fastq;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

public class FastqProcessor {

    // 假设conda安装路径，可根据实际情况调整
    private static final String CONDA_ENV = "common_env";
    private static final String HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz";
    private static final String HG38_FILE = "hg38.fa";
    private static final String FASTQ_DIR = "./fastq";
    private static final String BAM_DIR = "./bam";
    static boolean ignored;

    /**
     * 主流程：下载参考基因组、索引、处理所有fastq文件
     */
    public static void processFastqToBam() throws IOException, InterruptedException {
        // 确保输出目录存在
        ignored = new File(BAM_DIR).mkdirs();

        // 1. 准备参考基因组
        prepareReferenceGenome();

        // 2. 遍历fastq文件并处理
        File fastqFolder = new File(FASTQ_DIR);
        if (!fastqFolder.exists() || !fastqFolder.isDirectory()) {
            throw new IOException("Fastq目录不存在: " + FASTQ_DIR);
        }

        File[] fastqFiles = fastqFolder.listFiles((dir, name) ->
                name.endsWith(".fastq") || name.endsWith(".fq"));
        if (fastqFiles == null || fastqFiles.length == 0) {
            System.out.println("未找到fastq文件，跳过处理");
            return;
        }

        // 可选：识别配对端文件
        Map<String, List<File>> pairedFiles = detectPairedEnd(fastqFiles);

        for (Map.Entry<String, List<File>> entry : pairedFiles.entrySet()) {
            List<File> files = entry.getValue();
            String sampleName = entry.getKey();

            if (files.size() == 2) {
                // 双端测序
                processPairedEnd(files.get(0), files.get(1), sampleName);
            } else {
                // 单端测序
                processSingleEnd(files.get(0), sampleName);
            }
        }
    }

    /**
     * 下载并索引参考基因组
     */
    private static void prepareReferenceGenome() throws IOException, InterruptedException {
        File genomeFile = new File(HG38_FILE);
        if (!genomeFile.exists()) {
            System.out.println("下载hg38参考基因组...");
            // 下载并解压
            runCommand("wget", "-O", "hg38.fa.gz", HG38_URL);
            runCommand("gunzip", "hg38.fa.gz");
        } else {
            System.out.println("参考基因组已存在: " + HG38_FILE);
        }

        // 检查bwa索引文件
        File bwaIndex = new File(HG38_FILE + ".bwt");
        if (!bwaIndex.exists()) {
            System.out.println("构建bwa索引...");
            runCommandWithConda("bwa", "index", HG38_FILE);
        } else {
            System.out.println("bwa索引已存在");
        }
    }

    /**
     * 处理单端fastq文件
     */
    private static void processSingleEnd(File fastq, String sampleName) throws IOException, InterruptedException {
        String samFile = sampleName + ".sam";
        String unsortedBam = sampleName + ".bam";
        String sortedBam = sampleName + ".sorted.bam";
        File finalBam = new File(BAM_DIR, sortedBam);

        if (finalBam.exists()) {
            System.out.println("已存在排序bam文件，跳过: " + finalBam.getName());
            return;
        }

        System.out.println("处理单端样本: " + sampleName);
        // bwa mem 比对，通过 -o 直接输出 sam
        runCommandWithConda("bwa", "mem", "-t", "4", "-o", samFile, HG38_FILE, fastq.getAbsolutePath());
        // sam转bam
        runCommandWithConda("samtools", "view", "-bS", samFile, "-o", unsortedBam);
        // 排序bam
        runCommandWithConda("samtools", "sort", unsortedBam, "-o", finalBam.getAbsolutePath());
        // 清理中间文件
        ignored = new File(samFile).delete();
        ignored = new File(unsortedBam).delete();
        ignored = fastq.delete();
        System.out.println("完成: " + finalBam.getAbsolutePath());
    }

    /**
     * 处理双端fastq文件
     */
    private static void processPairedEnd(File fastq1, File fastq2, String sampleName) throws IOException, InterruptedException {
        String samFile = sampleName + ".sam";
        String unsortedBam = sampleName + ".bam";
        String sortedBam = sampleName + ".sorted.bam";
        File finalBam = new File(BAM_DIR, sortedBam);

        if (finalBam.exists()) {
            System.out.println("已存在排序bam文件，跳过: " + finalBam.getName());
            return;
        }

        System.out.println("处理双端样本: " + sampleName);
        runCommandWithConda("bwa", "mem", "-t", "4", "-o", samFile, HG38_FILE,
                fastq1.getAbsolutePath(), fastq2.getAbsolutePath());
        runCommandWithConda("samtools", "view", "-bS", samFile, "-o", unsortedBam);
        runCommandWithConda("samtools", "sort", unsortedBam, "-o", finalBam.getAbsolutePath());
        ignored = new File(samFile).delete();
        ignored = new File(unsortedBam).delete();
        ignored = fastq1.delete();
        ignored = fastq2.delete();
        System.out.println("完成: " + finalBam.getAbsolutePath());
    }

    /**
     * 检测配对端文件，返回样本名到文件列表的映射
     * 假设配对端文件命名模式为: 样本名_1.fastq 和 样本名_2.fastq
     */
    private static Map<String, List<File>> detectPairedEnd(File[] files) {
        Map<String, List<File>> map = new HashMap<>();
        Pattern p = Pattern.compile("^(.*)_([12])\\.(fastq|fq)$");
        for (File f : files) {
            String name = f.getName();
            java.util.regex.Matcher m = p.matcher(name);
            if (m.matches()) {
                String sample = m.group(1);
                map.computeIfAbsent(sample, k -> new ArrayList<>()).add(f);
            } else {
                // 单端文件，以文件名为样本名
                String sample = name.replaceFirst("\\.(fastq|fq)$", "");
                map.computeIfAbsent(sample, k -> new ArrayList<>()).add(f);
            }
        }
        return map;
    }

    /**
     * 执行系统命令（不使用conda环境）
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
     * 在conda环境中执行命令（使用conda run）
     */
    private static void runCommandWithConda(String... cmd) throws IOException, InterruptedException {
        // 构建命令: conda run -n ENV 实际命令...
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