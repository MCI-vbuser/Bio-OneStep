package com.vbuser;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class Main {

    private static final String CONDA_ENV_NAME = "common_env";
    private static final String PYTHON_VERSION = "3.13";
    private static final String[] PACKAGES = {
            "minimap", "fastp", "stringtie", "gffcompare", "suppa", "miranda", "gffread"
    };
    private static final String MINICONDA_INSTALL_PATH = System.getProperty("user.home") + "/miniconda3";
    private static final String MINICONDA_INSTALLER = "Miniconda3-latest-Linux-x86_64.sh";
    private static final String MINICONDA_URL = "https://repo.anaconda.com/miniconda/" + MINICONDA_INSTALLER;

    public static void main(String[] args) {
        envCheck();
    }

    private static void envCheck() {
        // 检查操作系统
        if (!System.getProperty("os.name").toLowerCase().contains("linux")) {
            System.out.println("警告：该程序主要针对 Linux 设计，在其他系统上可能运行异常。");
        }

        // 1. 获取 Conda 可执行文件路径
        String condaExec = findOrInstallConda();
        if (condaExec == null) {
            System.err.println("错误：无法找到或安装 Conda，请手动安装后重试。");
            return;
        }
        System.out.println("使用 Conda: " + condaExec);

        // 2. 初始化 Conda（新安装需要）
        initializeConda(condaExec);

        // 3. 确保 Conda 使用条款已被同意（显式执行 tos accept）
        if (!acceptCondaTOS(condaExec)) {
            System.err.println("错误：无法同意 Conda 使用条款，后续操作可能失败。");
            return;
        }

        // 4. 配置清华镜像（可选）
        runCommand(condaExec, "config", "--add", "channels", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/");

        // 5. 确保目标环境存在
        if (!condaEnvExists(condaExec)) {
            System.out.println("Conda 环境 " + CONDA_ENV_NAME + " 不存在，开始创建...");
            if (!createCondaEnv(condaExec)) {
                System.err.println("错误：创建 Conda 环境失败。");
                return;
            }
        } else {
            System.out.println("Conda 环境 " + CONDA_ENV_NAME + " 已存在。");
        }

        // 6. 安装依赖包
        if (!installPackages(condaExec)) {
            System.err.println("错误：安装依赖包失败。");
            return;
        }
        System.out.println("环境配置完成！");
    }

    /**
     * 初始化 Conda（执行 conda init）
     */
    private static void initializeConda(String condaExec) {
        System.out.println("初始化 Conda...");
        runShellCommandQuiet(condaExec + " init bash > /dev/null 2>&1");
    }

    /**
     * 查找系统 Conda，若找不到则静默安装 Miniconda
     */
    private static String findOrInstallConda() {
        if (isCondaAvailable("conda")) {
            return "conda";
        }

        String[] commonPaths = {
                MINICONDA_INSTALL_PATH + "/bin/conda",
                System.getProperty("user.home") + "/anaconda3/bin/conda",
                "/opt/anaconda3/bin/conda",
                "/usr/local/anaconda3/bin/conda"
        };
        for (String path : commonPaths) {
            if (new File(path).exists() && isCondaAvailable(path)) {
                return path;
            }
        }

        System.out.println("未找到 Conda，开始静默安装 Miniconda 到 " + MINICONDA_INSTALL_PATH);
        if (!installMiniconda()) {
            return null;
        }

        String installedConda = MINICONDA_INSTALL_PATH + "/bin/conda";
        return isCondaAvailable(installedConda) ? installedConda : null;
    }

    private static boolean isCondaAvailable(String condaCmd) {
        try {
            ProcessBuilder pb = new ProcessBuilder(condaCmd, "--version");
            pb.redirectErrorStream(true);
            Process p = pb.start();
            if (p.waitFor(10, TimeUnit.SECONDS) && p.exitValue() == 0) {
                return true;
            }
        } catch (Exception e) {
            // ignore
        }
        return false;
    }

    /**
     * 静默安装 Miniconda（无进度输出）
     */
    private static boolean installMiniconda() {
        String downloader;
        if (checkCommand("wget")) {
            downloader = "wget -q -O " + MINICONDA_INSTALLER + " " + MINICONDA_URL;
        } else if (checkCommand("curl")) {
            downloader = "curl -s -o " + MINICONDA_INSTALLER + " " + MINICONDA_URL;
        } else {
            System.err.println("错误：系统中未找到 wget 或 curl，无法下载 Miniconda 安装脚本。");
            return false;
        }

        System.out.println("下载 Miniconda 安装脚本...");
        if (!runShellCommandQuiet(downloader)) {
            System.err.println("下载失败，请检查网络连接。");
            return false;
        }
        System.out.println("下载完成。");

        System.out.println("运行安装脚本...");
        String installCmd = "bash " + MINICONDA_INSTALLER + " -b -p " + MINICONDA_INSTALL_PATH + " > /dev/null 2>&1";
        if (!runShellCommandQuiet(installCmd)) {
            System.err.println("安装失败，请检查磁盘空间和权限。");
            return false;
        }

        boolean ignored = new File(MINICONDA_INSTALLER).delete();
        createCondaConfig(); // 创建 .condarc 预先同意条款
        return true;
    }

    /**
     * 创建 .condarc 配置文件，预先同意条款
     */
    private static void createCondaConfig() {
        String condarcPath = System.getProperty("user.home") + "/.condarc";
        try (PrintWriter writer = new PrintWriter(new FileWriter(condarcPath))) {
            writer.println("channels:");
            writer.println("  - defaults");
            writer.println("anaconda_terms_of_service_accept: true");
            writer.println("auto_activate_base: false");
            System.out.println("已创建 .condarc 配置文件并预先同意条款。");
        } catch (IOException e) {
            System.err.println("警告：无法创建 .condarc 配置文件。");
        }
    }

    /**
     * 显式执行 conda tos accept，确保条款被接受
     */
    private static boolean acceptCondaTOS(String condaExec) {
        System.out.println("确保 Conda 使用条款已被同意...");

        // 尝试通过 conda config 设置（可能对旧版本有效）
        runShellCommandQuiet(condaExec + " config --set anaconda_terms_of_service_accept true > /dev/null 2>&1");

        // 显式执行 tos accept 命令，不带 --yes 参数
        String[] channels = {
                "https://repo.anaconda.com/pkgs/main",
                "https://repo.anaconda.com/pkgs/r"
        };
        boolean allAccepted = true;
        for (String channel : channels) {
            System.out.println("接受频道条款: " + channel);
            String cmd = condaExec + " tos accept --override-channels --channel " + channel;
            if (runShellCommand(cmd)) { // 使用 runShellCommand 以便看到输出
                System.err.println("警告：无法接受频道 " + channel + " 的条款");
                allAccepted = false;
            }
        }

        // 接受所有已配置的频道
        System.out.println("接受所有已配置频道的条款...");
        String cmdAll = condaExec + " tos accept";
        if (runShellCommand(cmdAll)) {
            System.err.println("警告：无法接受所有频道的条款");
            allAccepted = false;
        }

        return allAccepted;
    }

    /**
     * 执行 shell 命令，返回是否成功（不重定向输出，用于调试）
     */
    private static boolean runShellCommand(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", command);
            pb.redirectErrorStream(true);
            Process p = pb.start();
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println("[tos] " + line);
                }
            }
            int exitCode = p.waitFor();
            return exitCode != 0;
        } catch (Exception e) {
            System.err.println("执行命令失败: " + command);
            return true;
        }
    }

    /**
     * 静默执行 shell 命令
     */
    private static boolean runShellCommandQuiet(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", command);
            pb.redirectErrorStream(true);
            pb.redirectOutput(new File("/dev/null"));
            pb.redirectError(new File("/dev/null"));
            Process p = pb.start();
            int exitCode = p.waitFor();
            return exitCode == 0;
        } catch (Exception e) {
            return false;
        }
    }

    private static boolean checkCommand(String cmd) {
        try {
            ProcessBuilder pb = new ProcessBuilder("which", cmd);
            pb.redirectErrorStream(true);
            Process p = pb.start();
            return p.waitFor() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * 执行任意命令，设置环境变量，打印输出
     */
    private static boolean runCommand(String... cmdarray) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmdarray);
            pb.redirectErrorStream(true);
            pb.environment().put("CONDA_ALWAYS_YES", "true");
            pb.environment().put("CONDA_TERMS_OF_SERVICE_ACCEPT", "true");
            Process p = pb.start();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println("[exec] " + line);
                }
            }

            int exitCode = p.waitFor();
            if (exitCode != 0) {
                System.err.println("命令执行失败，退出码: " + exitCode + "，命令: " + String.join(" ", cmdarray));
                return false;
            }
            return true;
        } catch (Exception e) {
            System.err.println("执行命令时发生异常: " + String.join(" ", cmdarray));
            return false;
        }
    }

    /**
     * 静默执行命令，设置环境变量
     */
    private static boolean runCommandQuiet(String... cmdarray) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmdarray);
            pb.redirectOutput(new File("/dev/null"));
            pb.redirectError(new File("/dev/null"));
            pb.environment().put("CONDA_ALWAYS_YES", "true");
            pb.environment().put("CONDA_TERMS_OF_SERVICE_ACCEPT", "true");
            Process p = pb.start();
            int exitCode = p.waitFor();
            if (exitCode != 0) {
                System.err.println("命令执行失败，退出码: " + exitCode + "，命令: " + String.join(" ", cmdarray));
                return false;
            }
            return true;
        } catch (Exception e) {
            System.err.println("执行命令时发生异常: " + String.join(" ", cmdarray));
            return false;
        }
    }

    /**
     * 执行 conda 命令，自动添加 -y
     */
    private static boolean runCondaCommand(String condaExec, boolean quiet, String... args) {
        List<String> cmdList = new ArrayList<>();
        cmdList.add(condaExec);

        boolean yesAdded = false;
        for (String arg : args) {
            cmdList.add(arg);
            if ("-y".equals(arg) || "--yes".equals(arg)) {
                yesAdded = true;
            }
        }
        if (!yesAdded) {
            cmdList.add("-y");
        }

        if (quiet) {
            return runCommandQuiet(cmdList.toArray(new String[0]));
        } else {
            return runCommand(cmdList.toArray(new String[0]));
        }
    }

    private static boolean condaEnvExists(String condaExec) {
        try {
            ProcessBuilder pb = new ProcessBuilder(condaExec, "env", "list");
            pb.redirectErrorStream(true);
            pb.environment().put("CONDA_ALWAYS_YES", "true");
            pb.environment().put("CONDA_TERMS_OF_SERVICE_ACCEPT", "true");
            Process p = pb.start();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.trim().startsWith(CONDA_ENV_NAME + " ") || line.contains("/" + CONDA_ENV_NAME)) {
                        p.destroy();
                        return true;
                    }
                }
            }
            p.waitFor();
        } catch (Exception e) {
            System.out.println("检查 Conda 环境时发生异常: " + e.getMessage());
        }
        return false;
    }

    private static boolean createCondaEnv(String condaExec) {
        return runCondaCommand(condaExec, false, "create", "-n", CONDA_ENV_NAME, "python=" + PYTHON_VERSION);
    }

    private static boolean installPackages(String condaExec) {
        List<String> cmd = new ArrayList<>();
        cmd.add("install");
        cmd.add("-n");
        cmd.add(CONDA_ENV_NAME);
        cmd.add("-c");
        cmd.add("bioconda");
        Collections.addAll(cmd, PACKAGES);
        return runCondaCommand(condaExec, true, cmd.toArray(new String[0]));
    }
}