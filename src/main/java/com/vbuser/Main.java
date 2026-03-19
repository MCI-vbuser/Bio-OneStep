package com.vbuser;

import java.io.*;
import java.util.*;
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

    // 用于存放源码安装的 sra-toolkit bin 目录，便于后续命令直接使用
    private static String extraBinPath = null;

    static boolean ignored;

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

        // 7. 检查并安装 sra-toolkit
        if (!checkAndInstallSraToolkit()) {
            System.err.println("警告：sra-toolkit 安装失败，后续可能需要手动安装。");
        } else {
            System.out.println("sra-toolkit 已就绪。");
        }

        System.out.println("环境配置完成！");
    }

    // ==================== 新增：sra-toolkit 检测与安装 ====================

    /**
     * 检查 sra-toolkit 是否可用，若不可用则尝试安装（apt 或源码）
     */
    private static boolean checkAndInstallSraToolkit() {
        String userHome = System.getProperty("user.home");

        // 检查是否已有通过源码安装的版本（例如上次运行留下的 sratoolkit-current）
        File possibleBin = new File(userHome + "/sratoolkit-current/bin/prefetch");
        if (possibleBin.exists() && possibleBin.canExecute()) {
            extraBinPath = userHome + "/sratoolkit-current/bin";
            System.out.println("sra-toolkit 已安装 (通过源码，位于 " + extraBinPath + ")");
            return true;
        }

        // 检查系统 PATH 中是否有 prefetch
        if (checkCommand("prefetch")) {
            System.out.println("sra-toolkit 已安装 (prefetch 在系统 PATH 中)。");
            return true;
        }

        System.out.println("未找到 sra-toolkit，开始安装...");

        // 尝试使用 sudo apt 安装
        if (attemptSudoAptInstall()) {
            // 再次检查 prefetch
            if (checkCommand("prefetch")) {
                System.out.println("通过 apt 安装 sra-toolkit 成功。");
                return true;
            } else {
                System.out.println("apt 安装可能未成功，将尝试源码安装。");
            }
        } else {
            System.out.println("无法使用 apt 安装，将尝试源码安装。");
        }

        // 源码安装
        return installSraToolkitFromSource();
    }

    /**
     * 尝试使用 sudo apt 安装 sra-toolkit（支持交互式密码输入）
     */
    private static boolean attemptSudoAptInstall() {
        boolean isRoot = "root".equals(System.getProperty("user.name"));
        if (isRoot) {
            System.out.println("当前用户为 root，直接执行 apt install...");
            return runCommand("apt", "install", "-y", "sra-toolkit");
        } else {
            // 先尝试无密码 sudo
            if (runShellCommandQuiet("sudo -n true")) {
                System.out.println("检测到无密码 sudo，尝试使用 sudo apt install...");
                return runCommand("sudo", "apt", "install", "-y", "sra-toolkit");
            } else {
                // 需要密码：使用 inheritIO 让用户直接在终端输入密码
                System.out.println("需要密码的 sudo，正在弹出密码提示（请在终端中输入密码）...");
                try {
                    ProcessBuilder pb = new ProcessBuilder("sudo", "apt", "install", "-y", "sra-toolkit");
                    pb.inheritIO(); // 将子进程的输入输出连接到当前 Java 进程的终端
                    Process p = pb.start();
                    int exitCode = p.waitFor();
                    if (exitCode == 0) {
                        return true;
                    } else {
                        System.err.println("sudo apt 安装失败，退出码: " + exitCode);
                        return false;
                    }
                } catch (Exception e) {
                    System.err.println("执行 sudo apt 时发生异常: " + e.getMessage());
                    return false;
                }
            }
        }
    }

    /**
     * 从源码编译安装 sra-toolkit：下载、解压、创建软链接，并设置 extraBinPath
     */
    private static boolean installSraToolkitFromSource() {
        String userHome = System.getProperty("user.home");
        String downloadUrl = "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz";
        String tarFileName = "sratoolkit.current-ubuntu64.tar.gz";

        // 检查下载工具
        boolean useWget = checkCommand("wget");
        boolean useCurl = checkCommand("curl");
        if (!useWget && !useCurl) {
            System.err.println("错误：系统中未找到 wget 或 curl，无法下载 sra-toolkit。");
            return false;
        }

        String downloadCmd;
        if (useWget) {
            downloadCmd = "wget -q -O " + tarFileName + " " + downloadUrl;
        } else {
            downloadCmd = "curl -s -o " + tarFileName + " " + downloadUrl;
        }

        System.out.println("下载 sra-toolkit 源码包...");
        if (!runShellCommandQuiet(downloadCmd)) {
            System.err.println("下载失败，请检查网络连接。");
            return false;
        }
        System.out.println("下载完成。");

        // 解压到用户主目录
        System.out.println("解压中...");
        String extractCmd = "tar -xzf " + tarFileName + " -C " + userHome;
        if (!runShellCommandQuiet(extractCmd)) {
            System.err.println("解压失败。");
            ignored = new File(tarFileName).delete();
            return false;
        }
        System.out.println("解压完成。");

        // 查找解压后的目录名（格式如 sratoolkit.3.0.0-ubuntu64）
        File homeDir = new File(userHome);
        File[] sraDirs = homeDir.listFiles((dir, name) -> name.startsWith("sratoolkit.") && name.endsWith("-ubuntu64"));
        if (sraDirs == null || sraDirs.length == 0) {
            System.err.println("无法找到解压后的 sratoolkit 目录。");
            ignored = new File(tarFileName).delete();
            return false;
        }
        File sraDir = sraDirs[0]; // 通常只有一个
        String sraDirPath = sraDir.getAbsolutePath();

        // 创建软链接 ~/sratoolkit-current 指向该目录
        String linkPath = userHome + "/sratoolkit-current";
        File linkFile = new File(linkPath);
        if (linkFile.exists()) {
            if (!linkFile.delete()) {
                System.err.println("警告：无法删除已有的软链接 " + linkPath);
            }
        }
        String lnCmd = "ln -s " + sraDirPath + " " + linkPath;
        if (!runShellCommandQuiet(lnCmd)) {
            System.err.println("创建软链接失败。");
            ignored = new File(tarFileName).delete();
            return false;
        }

        // 设置 extraBinPath 为该软链接下的 bin 目录
        extraBinPath = linkPath + "/bin";

        // 清理下载的 tar 文件
        ignored = new File(tarFileName).delete();

        // 验证 prefetch 是否可用（此时 extraBinPath 已生效）
        System.out.println("验证安装...");
        if (runCommand("prefetch", "--version")) {
            System.out.println("sra-toolkit 源码安装成功。");
            return true;
        } else {
            System.err.println("警告：prefetch 命令仍不可用，安装可能有问题。");
            return false;
        }
    }

    // ==================== 修改原有命令执行方法，支持 extraBinPath ====================

    /**
     * 执行任意命令，设置环境变量，打印输出（已支持 extraBinPath）
     */
    private static boolean runCommand(String... cmdarray) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmdarray);
            // 将 extraBinPath 添加到 PATH 环境变量
            if (extraBinPath != null) {
                Map<String, String> env = pb.environment();
                env.merge("PATH", extraBinPath, (a, b) -> b + ":" + a);
            }
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
     * 静默执行命令，设置环境变量（已支持 extraBinPath）
     */
    private static boolean runCommandQuiet(String... cmdarray) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmdarray);
            if (extraBinPath != null) {
                Map<String, String> env = pb.environment();
                env.merge("PATH", extraBinPath, (a, b) -> b + ":" + a);
            }
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
     * 执行 shell 命令，返回是否成功（不重定向输出，用于调试）（已支持 extraBinPath）
     */
    private static boolean runShellCommand(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", command);
            if (extraBinPath != null) {
                Map<String, String> env = pb.environment();
                env.merge("PATH", extraBinPath, (a, b) -> b + ":" + a);
            }
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
     * 静默执行 shell 命令（已支持 extraBinPath）
     */
    private static boolean runShellCommandQuiet(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", command);
            if (extraBinPath != null) {
                Map<String, String> env = pb.environment();
                env.merge("PATH", extraBinPath, (a, b) -> b + ":" + a);
            }
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

    // ==================== 以下为原有方法（未修改，仅保持完整性） ====================

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