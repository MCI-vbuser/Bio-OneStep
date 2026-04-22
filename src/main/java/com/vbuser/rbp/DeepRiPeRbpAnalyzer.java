package com.vbuser.rbp;

import java.io.*;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DeepRiPeRbpAnalyzer {

    private static final String CONDA_ENV = "deepripe_env";
    private static final String SCRIPT_DIR = "./scripts/rbp";
    private static final String MODEL_DIR = SCRIPT_DIR;
    private static final String INPUT_FASTA = "./as/novel.fa";
    private static final String OUTPUT_BASE_DIR = "./rbp";

    private static final String OUTPUT_CSV = OUTPUT_BASE_DIR + "/rbp_predictions.csv";
    private static final String FILTERED_CSV = OUTPUT_BASE_DIR + "/rbp_predictions_binding_true.csv";
    private static final String REPORTS_DIR = OUTPUT_BASE_DIR + "/reports";
    private static final String ADVANCED_DIR = OUTPUT_BASE_DIR + "/advanced_analysis";

    public static void runFullAnalysis() throws IOException, InterruptedException {
        System.out.println("=== DeepRiPe RBP 分析开始 ===");
        Files.createDirectories(Paths.get(OUTPUT_BASE_DIR));
        Files.createDirectories(Paths.get(REPORTS_DIR));
        Files.createDirectories(Paths.get(ADVANCED_DIR));
        // 确保 motif 目录存在（用于存放报告）
        Files.createDirectories(Paths.get("./rbp"));
        if(new File("./rbp/rbp_data.db").exists()){
            System.out.println("RBP分析已完成,跳过...");
            return;
        }

        setupEnvironment();
        runPrediction();
        runFilter();
        runVisualization();
        runAdvancedAnalysis();
        runRbpReport();          // 新增：生成交互式 HTML 报告
        System.out.println("=== DeepRiPe RBP 分析完成，结果保存在 " + OUTPUT_BASE_DIR + " ===");
    }

    private static void setupEnvironment() throws IOException, InterruptedException {
        // 检查环境是否存在
        ProcessBuilder checkPb = new ProcessBuilder("conda", "env", "list");
        checkPb.redirectErrorStream(true);
        Process checkProcess = checkPb.start();
        String output = readStream(checkProcess.getInputStream());
        boolean exists = output.contains(CONDA_ENV);
        checkProcess.waitFor();

        if (!exists) {
            System.out.println("创建 Conda 环境: " + CONDA_ENV);
            runCommand("conda", "create", "-n", CONDA_ENV, "python=3.8", "-y");
        } else {
            System.out.println("Conda 环境已存在: " + CONDA_ENV);
        }

        // 安装 Python 依赖（全部使用 pip，避免 conda 依赖冲突）
        System.out.println("安装 Python 依赖...");
        String[] deps = {
                "tensorflow==2.13.0", "pandas", "numpy", "h5py",
                "matplotlib", "seaborn", "biopython", "networkx",
                "plotly", "scipy"          // 新增 plotly 和 scipy（用于 KDE）
        };
        for (String dep : deps) {
            runCondaCommand("pip", "install", dep);
        }

        System.out.println("环境准备完成");
    }

    private static void runPrediction() throws IOException, InterruptedException {
        System.out.println("运行 RBP 结合预测...");
        File inputFile = new File(INPUT_FASTA);
        if (!inputFile.exists()) {
            throw new IOException("输入文件不存在: " + INPUT_FASTA);
        }
        ensureModelFiles();

        String scriptPath = Paths.get(SCRIPT_DIR, "rbp_inference.py").toString();
        runCondaCommand(
                "python", scriptPath,
                "--input", INPUT_FASTA,
                "--input-type", "fasta",
                "--output", OUTPUT_CSV,
                "--model-dir", MODEL_DIR
        );
        System.out.println("预测完成，结果保存至: " + OUTPUT_CSV);
    }

    private static void runFilter() throws IOException, InterruptedException {
        System.out.println("过滤预测结果...");
        String scriptPath = Paths.get(SCRIPT_DIR, "filter.py").toString();
        // 注意：filter.py 的输入是位置参数，不是 --input
        runCondaCommand(
                "python", scriptPath,
                OUTPUT_CSV,                     // 位置参数
                "--output", FILTERED_CSV,
                "--min-score", "0.5"
        );
        System.out.println("过滤完成，结果保存至: " + FILTERED_CSV);
    }

    private static void runVisualization() throws IOException, InterruptedException {
        System.out.println("生成可视化图表...");
        String scriptPath = Paths.get(SCRIPT_DIR, "visualize_results.py").toString();
        runCondaCommand(
                "python", scriptPath,
                "--input", FILTERED_CSV,
                "--output", REPORTS_DIR,
                "--create-report"
        );
        System.out.println("可视化完成，报告保存至: " + REPORTS_DIR);
    }

    private static void runAdvancedAnalysis() throws IOException, InterruptedException {
        System.out.println("运行高级分析...");
        String scriptPath = Paths.get(SCRIPT_DIR, "advanced_analysis.py").toString();
        runCondaCommand(
                "python", scriptPath,
                "--results", FILTERED_CSV,
                "--output", ADVANCED_DIR,
                "--network", "--co-binding"
        );
        System.out.println("高级分析完成，结果保存至: " + ADVANCED_DIR);
    }

    /**
     * 新增：生成交互式 RBP 报告（Plotly）
     */
    private static void runRbpReport() throws IOException, InterruptedException {
        System.out.println("生成 RBP 交互式报告...");
        File filteredCsv = new File(FILTERED_CSV);
        if (!filteredCsv.exists()) {
            System.err.println("警告：过滤后的 CSV 不存在，跳过报告生成。");
            return;
        }
        String scriptPath = Paths.get(SCRIPT_DIR, "rbp.py").toString();
        File scriptFile = new File(scriptPath);
        if (!scriptFile.exists()) {
            System.err.println("警告：rbp.py 脚本不存在，跳过报告生成。");
            return;
        }
        String reportHtml = "./rbp/rbp_report.html";
        String dbPath = "./rbp/rbp_data.db";
        runCondaCommand(
                "python", scriptPath,
                "--input", FILTERED_CSV,
                "--output", reportHtml,
                "--db", dbPath
        );
        System.out.println("RBP 报告已生成: " + reportHtml);
        System.out.println("RBP 数据库已生成: " + dbPath);
    }

    private static void ensureModelFiles() throws IOException {
        String[] models = {"model_RBPshigh.h5", "model_RBPsmed.h5", "model_RBPslow.h5"};
        for (String model : models) {
            Path modelPath = Paths.get(MODEL_DIR, model);
            if (!Files.exists(modelPath)) {
                throw new IOException("模型文件缺失: " + modelPath.toAbsolutePath() +
                        "\n请将预训练权重文件放入 " + MODEL_DIR + " 目录");
            }
        }
    }

    private static void runCommand(String... args) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder(args);
        pb.redirectErrorStream(true);
        System.out.println("执行命令: " + String.join(" ", args));
        Process process = pb.start();
        String output = readStream(process.getInputStream());
        int exitCode = process.waitFor();
        if (exitCode != 0) {
            throw new IOException("命令执行失败，退出码 " + exitCode + "\n" + output);
        }
    }

    private static void runCondaCommand(String... args) throws IOException, InterruptedException {
        List<String> fullCmd = new ArrayList<>();
        fullCmd.add("conda");
        fullCmd.add("run");
        fullCmd.add("-n");
        fullCmd.add(DeepRiPeRbpAnalyzer.CONDA_ENV);
        fullCmd.addAll(Arrays.asList(args));
        runCommand(fullCmd.toArray(new String[0]));
    }

    private static String readStream(InputStream is) throws IOException {
        StringBuilder sb = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
            String line;
            while ((line = reader.readLine()) != null) {
                sb.append(line).append("\n");
            }
        }
        return sb.toString();
    }
}