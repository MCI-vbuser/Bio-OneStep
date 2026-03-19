package com.vbuser.pre;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import static com.vbuser.Main.findOrInstallConda;
import static com.vbuser.Main.runCondaCommand;

public class Meme {
    /**
     * 检查指定名称的 Conda 环境是否存在
     */
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
                    if (line.trim().startsWith("meme_env" + " ") || line.contains("/" + "meme_env")) {
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

    /**
     * 创建指定名称和 Python 版本的 Conda 环境
     */
    private static boolean createCondaEnv(String condaExec) {
        return runCondaCommand(condaExec, false, "create", "-n", "meme_env", "python=" + "3.6");
    }

    /**
     * 检查指定 Conda 环境中是否已安装某个包
     */
    private static boolean isPackageInstalled(String condaExec) {
        return runCondaCommand(condaExec, true, "list", "-n", "meme_env", "meme");
    }

    /**
     * 在 meme_env 环境中安装 meme
     */
    private static boolean installMeme(String condaExec) {
        System.out.println("安装 meme 到 meme_env...");
        return runCondaCommand(condaExec, true, "install", "-n", "meme_env", "-c", "bioconda", "meme");
    }

    /**
     * 创建 meme_env 环境（Python 3.6）并安装 meme
     * 若环境已存在则检查 meme 是否已安装，若未安装则进行安装
     * @return true 表示环境准备成功，false 表示失败
     */
    public static boolean createMemeEnv() {
        String condaExec = findOrInstallConda();
        if (condaExec == null) {
            System.err.println("错误：无法找到或安装 Conda，无法创建 meme_env。");
            return false;
        }

        if (condaEnvExists(condaExec)) {
            System.out.println("环境 meme_env 已存在。");
            if (!isPackageInstalled(condaExec)) {
                System.out.println("meme 未安装，开始安装...");
                return installMeme(condaExec);
            } else {
                System.out.println("meme 已安装，无需操作。");
                return true;
            }
        } else {
            System.out.println("创建环境 meme_env (Python 3.6)...");
            if (!createCondaEnv(condaExec)) {
                System.err.println("错误：创建环境失败。");
                return false;
            }
            return installMeme(condaExec);
        }
    }
}
