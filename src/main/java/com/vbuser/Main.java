package com.vbuser;

import com.vbuser.gtf.SingleGTF;
import com.vbuser.pre.CommonEnv;
import com.vbuser.pre.DownloadSRA;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {

    private static final String CONDA_ENV = "common_env";   // 与 SingleGTF 保持一致

    public static void main(String[] args) {
        try {
            WebConsole.start();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        CommonEnv.envCheck();
        cls();
        System.out.println("环境检查通过");

        // 确保参考基因组已准备就绪
        try {
            ensureReferenceGenome();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }

        prepareDataAndRunGTF();

        System.out.println("战至最后一刻!自刎归天!");
        try {
            Thread.sleep(1000);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
        WebConsole.stopServer();
        WebConsole.killProcess();
    }

    /**
     * 确保 hg38.fa 参考基因组及其 minimap2 索引存在
     */
    private static void ensureReferenceGenome() throws IOException, InterruptedException {
        File genomeFile = new File("hg38.fa");
        if (!genomeFile.exists()) {
            System.out.println("下载 hg38 参考基因组...");
            // 使用系统命令 wget 和 gunzip
            runSystemCommand("wget", "-O", "hg38.fa.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz");
            runSystemCommand("gunzip", "hg38.fa.gz");
        } else {
            System.out.println("参考基因组已存在: hg38.fa");
        }

        // 检查 minimap2 索引
        File mmiIndex = new File("hg38.fa.mmi");
        if (!mmiIndex.exists()) {
            System.out.println("构建 minimap2 索引...");
            runCommandWithConda();
        } else {
            System.out.println("minimap2 索引已存在: hg38.fa.mmi");
        }
    }

    private static void prepareDataAndRunGTF() {
        System.out.println("[1] 使用本地 SRA 文件");
        System.out.println("[2] 下载新的 SRA 文件");
        System.out.println("[3] 已有 BAM 文件（跳过 SRA 转换和比对）");
        String ans = WebConsole.readLine();

        File input;
        boolean sortedEnsured = false;

        switch (ans) {
            case "1": {
                System.out.println("请指定 SRA 文件或目录路径：");
                String path = WebConsole.readLine();
                input = new File(path);
                break;
            }
            case "2":
                System.out.println("请指定 SRR 编号：");
                String srr = WebConsole.readLine();
                DownloadSRA.download(srr);
                input = new File("./sra");
                break;
            case "3": {
                System.out.println("请指定 BAM 文件或目录路径：");
                String path = WebConsole.readLine();
                if (path.contains(" -y")) {
                    sortedEnsured = true;
                    path = path.split(" -y")[0];
                }
                input = new File(path);
                break;
            }
            default:
                System.out.println("无效选择，退出");
                return;
        }

        // 调用 SingleGTF 处理
        SingleGTF.handle(input, sortedEnsured);
    }

    private static void cls() {
        try {
            new ProcessBuilder("bash", "-c", "clear").inheritIO().start().waitFor();
        } catch (InterruptedException | IOException e) {
            throw new RuntimeException(e);
        }
    }

    // ---------------------- 命令执行辅助方法 ----------------------
    private static void runSystemCommand(String... cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        Process process = pb.start();
        String output = readStream(process.getInputStream());
        int exitCode = process.waitFor();
        if (exitCode != 0) {
            throw new IOException("命令执行失败: " + String.join(" ", cmd) + "\n" + output);
        }
    }

    private static void runCommandWithConda() throws IOException, InterruptedException {
        List<String> fullCmd = new ArrayList<>();
        fullCmd.add("conda");
        fullCmd.add("run");
        fullCmd.add("-n");
        fullCmd.add(CONDA_ENV);
        fullCmd.addAll(Arrays.asList("minimap2", "-d", "hg38.fa.mmi", "hg38.fa"));
        runSystemCommand(fullCmd.toArray(new String[0]));
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