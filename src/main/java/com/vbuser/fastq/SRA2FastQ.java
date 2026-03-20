package com.vbuser.fastq;

import java.io.File;

public class SRA2FastQ {

    static boolean ignored;

    public static void handle(File rootDir) {
        // 递归遍历所有文件和目录
        processSraFiles(rootDir);
    }

    /**
     * 递归处理目录下的所有 .sra 文件
     * @param file 文件或目录
     */
    private static void processSraFiles(File file) {
        if (file == null) return;

        if (file.isFile() && file.getName().endsWith(".sra")) {
            // 处理单个 .sra 文件
            String cmd = "fastq-dump --split-files " + file.getAbsolutePath() + " -O ./fastq";
            // 创建输出目录（只需创建一次，但多次创建也无妨）
            ignored = new File("./fastq").mkdirs();
            try {
                Process process = Runtime.getRuntime().exec(cmd);
                process.waitFor();
                // 删除原文件（可选，根据需求决定）
                ignored = file.delete();
            } catch (Exception e) {
                throw new RuntimeException("处理文件失败: " + file.getAbsolutePath(), e);
            }
        } else if (file.isDirectory()) {
            // 递归遍历子目录
            File[] children = file.listFiles();
            if (children != null) {
                for (File child : children) {
                    processSraFiles(child);
                }
            }
        }
    }

}
