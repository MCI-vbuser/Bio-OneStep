package com.vbuser;

import com.vbuser.apa.APAtrapPipeline;
import com.vbuser.gtf.SingleGTF;
import com.vbuser.pre.CommonEnv;
import com.vbuser.pre.DownloadSRA;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
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

        // ========================================
        System.out.println("开始 AS 事件分析...");
        try {
            runAsPipeline();
        } catch (IOException | InterruptedException e) {
            System.err.println("AS 分析失败: " + e.getMessage());
            // 可根据需要决定是否继续
        }
        // ========================================
        System.out.println("开始 APA 分析...");
        try {
            APAtrapPipeline.run(new File("."));
        } catch (Exception e) {
            System.err.println("APA 分析失败: " + e.getMessage());
        }
        // ========================================
        try {
            runMotifAnalysis();
        } catch (Exception e) {
            System.err.println("Motif 分析失败: " + e.getMessage());
        }
        // ========================================
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
            runSystemCommand("wget", "-O", "hg38.fa.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz");
            runSystemCommand("gunzip", "hg38.fa.gz");
        } else {
            System.out.println("参考基因组已存在: hg38.fa");
        }

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

    private static void runAsPipeline() throws IOException, InterruptedException {
        String baseDir = new File(".").getAbsolutePath();
        String outputDir = new File("./as").getAbsolutePath();
        String scriptPath = "scripts/as/run_as_pipeline.py";

        if(new File(outputDir,"stats.db").exists()) {
            System.out.println("AS已完成，跳过分析");
            return;
        }

        File scriptFile = new File(scriptPath);
        if (!scriptFile.exists()) {
            throw new IOException("Python 脚本不存在: " + scriptFile.getAbsolutePath());
        }

        List<String> cmd = new ArrayList<>();
        cmd.add("conda");
        cmd.add("run");
        cmd.add("--no-capture-output");
        cmd.add("-n");
        cmd.add(CONDA_ENV);
        cmd.add("python");
        cmd.add(scriptPath);
        cmd.add("--base_dir");
        cmd.add(baseDir);
        cmd.add("--output_dir");
        cmd.add(outputDir);

        System.out.println("执行命令: " + String.join(" ", cmd));
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.environment().put("PYTHONUNBUFFERED", "1");  // 强制 Python 无缓冲输出
        pb.redirectErrorStream(true);
        Process process = pb.start();

        Thread outputReader = new Thread(() -> {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
            } catch (IOException e) {
                System.err.println("读取子进程输出时出错: " + e.getMessage());
            }
        });
        outputReader.setDaemon(true);
        outputReader.start();

        int exitCode = process.waitFor();
        outputReader.join(1000);
        if (exitCode != 0) {
            throw new IOException("AS 分析脚本执行失败，退出码: " + exitCode);
        }
        System.out.println("AS 分析完成。");
    }

    private static void runMotifAnalysis() throws IOException, InterruptedException {
        System.out.println("开始 Motif 分析...");

        if(new File("./motif/motifs.db").exists()) {
            System.out.println("Motif分析已完成,跳过...");
            return;
        }

        // ---------- 1. 确保参考基序数据库（RNA → DNA 转换） ----------
        File rbpDbRna = new File("Ray2013_rbp_RNA.meme");
        if (!rbpDbRna.exists()) {
            System.out.println("下载 Ray2013 RBP motif 数据库（RNA 格式）...");
            runSystemCommand("wget", "-O", "Ray2013_rbp_RNA.meme",
                    "https://raw.githubusercontent.com/xypan1232/iDeep/master/Ray2013_rbp_RNA.meme");
        }

        File rbpDbDna = new File("Ray2013_rbp_RNA.meme.dna");
        if (!rbpDbDna.exists()) {
            System.out.println("将 RNA motif 数据库转换为 DNA 格式（替换 U -> T）...");
            try {
                // 尝试用 sed（若系统支持且 ProcessBuilder 可重定向，通常失败回退 Java）
                runSystemCommand("sed", "s/U/T/g", rbpDbRna.getAbsolutePath(), ">", rbpDbDna.getAbsolutePath());
            } catch (IOException e) {
                System.out.println("sed 重定向失败，改用 Java 文本替换...");
                convertRnaMotifToDna(rbpDbRna, rbpDbDna);
            }
        }

        // ---------- 2. 确保 novel.fa（MSTRG 转录本序列）存在 ----------
        File novelFa = ensureNovelFasta();

        // ---------- 3. 背景模型 ----------
        File backgroundModel = new File("./apa/background.model");
        if (!backgroundModel.exists()) {
            System.out.println("背景模型不存在，从输入序列生成...");
            generateBackgroundModel(novelFa, backgroundModel);
        }

        // ---------- 4. 创建输出目录 ----------
        File motifDir = new File("./motif");
        if (!motifDir.exists() && !motifDir.mkdirs()) {
            throw new IOException("无法创建目录: " + motifDir);
        }

        String condaMeme = "meme_env";

        // ---------- 5. 运行 MEME ----------
        System.out.println("运行 MEME...");
        List<String> memeArgs = new ArrayList<>(Arrays.asList(
                "meme", novelFa.getAbsolutePath(),
                "-oc", new File(motifDir, "meme_results").getAbsolutePath(),
                "-dna", "-mod", "anr", "-nmotifs", "10", "-minw", "6", "-maxw", "20"
        ));
        if (backgroundModel.exists()) {
            memeArgs.add("-bfile");
            memeArgs.add(backgroundModel.getAbsolutePath());
        }
        runCondaCommand(condaMeme, memeArgs.toArray(new String[0]));

        // ---------- 6. 运行 DREME ----------
        System.out.println("运行 DREME...");
        runCondaCommand(condaMeme,
                "dreme-py3", "-p", novelFa.getAbsolutePath(),
                "-oc", new File(motifDir, "dreme_results").getAbsolutePath(),
                "-dna", "-e", "0.05", "-m", "5");

        // ---------- 7. 运行 Tomtom ----------
        System.out.println("运行 Tomtom...");
        File memeTxt = new File(motifDir, "meme_results/meme.txt");
        File dremeTxt = new File(motifDir, "dreme_results/dreme.txt");
        List<String> tomtomInputs = new ArrayList<>();
        if (memeTxt.exists()) tomtomInputs.add(memeTxt.getAbsolutePath());
        if (dremeTxt.exists()) tomtomInputs.add(dremeTxt.getAbsolutePath());

        if (!tomtomInputs.isEmpty()) {
            List<String> tomtomCmd = new ArrayList<>(Arrays.asList(
                    "tomtom", "-oc", new File(motifDir, "tomtom_results").getAbsolutePath()
            ));
            tomtomCmd.addAll(tomtomInputs);
            tomtomCmd.add(rbpDbDna.getAbsolutePath());
            runCondaCommand(condaMeme, tomtomCmd.toArray(new String[0]));
        } else {
            System.err.println("警告：未找到 meme.txt 或 dreme.txt，跳过 Tomtom。");
        }

        // ---------- 8. 生成报告 ----------
        System.out.println("生成 Motif 统计报告...");
        File motifPy = new File("./scripts/motif/motif.py");
        if (!motifPy.exists()) {
            System.err.println("警告：motif.py 脚本缺失，跳过报告生成。");
            return;
        }



        System.out.println("Motif 分析完成，结果保存于 " + motifDir.getAbsolutePath());
    }

    /**
     * 确保 novel.fa 文件存在，若不存在则通过 gffread 从增强 GTF 生成
     */
    private static File ensureNovelFasta() throws IOException, InterruptedException {
        // 尝试多个可能的路径
        File novelFa = new File("./as/novel.fa");
        if (novelFa.exists()) return novelFa;

        novelFa = new File("./as/novel_sequences.fa");
        if (novelFa.exists()) return novelFa;

        novelFa = new File("./as/novel_transcripts.fa");
        if (novelFa.exists()) return novelFa;

        // 都不存在，则生成
        System.out.println("未找到 MSTRG 序列文件，现在从增强 GTF 生成 novel.fa ...");
        File mergedGtf = new File("./as/gencode.v45.enhanced.gtf");
        if (!mergedGtf.exists()) {
            throw new IOException("增强 GTF 文件不存在，无法生成 novel.fa: " + mergedGtf.getAbsolutePath());
        }
        File refGenome = new File("hg38.fa");
        if (!refGenome.exists()) {
            throw new IOException("参考基因组 hg38.fa 不存在，无法提取序列。");
        }

        // 输出文件路径
        novelFa = new File("./as/novel.fa");

        // 使用 gffread 提取所有转录本 cDNA 序列（通过 conda 环境）
        System.out.println("运行 gffread 提取转录本序列...");
        runCondaCommand(CONDA_ENV, "gffread", "-w", novelFa.getAbsolutePath(),
                "-g", refGenome.getAbsolutePath(), mergedGtf.getAbsolutePath());

        // 可选：过滤只保留 MSTRG 开头的序列
        File filteredFa = new File("./as/novel_mstrg_only.fa");
        filterFastaByHeaderPrefix(novelFa, filteredFa);

        if (filteredFa.length() > 0) {
            Files.move(filteredFa.toPath(), novelFa.toPath(), StandardCopyOption.REPLACE_EXISTING);
        }
        System.out.println("novel.fa 生成完成: " + novelFa.getAbsolutePath());
        return novelFa;
    }

    /**
     * 过滤 FASTA 文件，仅保留 header 以指定前缀开头的序列
     */
    private static void filterFastaByHeaderPrefix(File input, File output) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(input));
             BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
            String line;
            boolean keep = false;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    keep = line.substring(1).trim().startsWith("MSTRG");
                }
                if (keep) {
                    writer.write(line + "\n");
                }
            }
        }
    }

    /**
     * 将 RNA motif 数据库（含字母 U）转换为 DNA 版本（U -> T）
     */
    private static void convertRnaMotifToDna(File input, File output) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(input));
             BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
            String line;
            while ((line = reader.readLine()) != null) {
                writer.write(line.replace('U', 'T').replace('u', 't'));
                writer.newLine();
            }
        }
    }

    /**
     * 使用 MEME suite 的 fasta-get-markov 生成背景模型
     */
    private static void generateBackgroundModel(File inputFasta, File outputModel) throws IOException, InterruptedException {
        runCondaCommand("meme_env", "fasta-get-markov", "-m", "0", inputFasta.getAbsolutePath(), outputModel.getAbsolutePath());
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

    private static void runCondaCommand(String env, String... args) throws IOException, InterruptedException {
        List<String> cmd = new ArrayList<>();
        cmd.add("conda");
        cmd.add("run");
        cmd.add("-n");
        cmd.add(env);
        cmd.addAll(Arrays.asList(args));
        runSystemCommand(cmd.toArray(new String[0]));
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