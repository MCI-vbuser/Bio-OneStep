package com.vbuser.pre;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DownloadSRA {
    public static void download(String sra_txt) {
        String[] sra = phrase(sra_txt.replaceAll(" ", ""));
        System.out.println("是否确认下载如下SRA文件: " + Arrays.toString(sra) + "? [y|n]");
        String input = System.console().readLine();
        if ("y".equals(input)) {
            String sraString = String.join("\n", sra);
            writeString2Txt(sraString);
            try {
                // 使用 ProcessBuilder 更安全地处理命令和参数
                ProcessBuilder pb = new ProcessBuilder(
                        "prefetch","--progress", "--option-file", "./sra.txt", "-O", "./sra"
                );
                pb.inheritIO(); // 将子进程的输入输出重定向到当前进程，便于观察
                Process process = pb.start();
                int exitCode = process.waitFor();
                if (exitCode != 0) {
                    System.err.println("prefetch 命令执行失败，退出码: " + exitCode);
                }
                // 无论成功与否，尝试删除临时文件
                File tempFile = new File("./sra.txt");
                if (tempFile.exists() && !tempFile.delete()) {
                    System.err.println("警告：无法删除临时文件 ./sra.txt");
                }
            } catch (IOException | InterruptedException e) {
                throw new RuntimeException("执行 prefetch 命令时出错", e);
            }
        }
    }

    private static String[] phrase(String sra) {
        List<String> result = new ArrayList<>();
        String[] parts = sra.split(",");
        for (String part : parts) {
            if (part.contains("-")) {
                // 处理范围，例如 SRR123-125
                String[] range = part.split("-");
                if (range.length != 2) {
                    // 格式错误，忽略或抛出异常，这里直接跳过
                    continue;
                }
                String startStr = range[0].replace("SRR", "");
                String endStr = range[1].replace("SRR", "");
                try {
                    int start = Integer.parseInt(startStr);
                    int end = Integer.parseInt(endStr);
                    for (int i = start; i <= end; i++) {
                        result.add("SRR" + i);
                    }
                } catch (NumberFormatException e) {
                    // 解析数字失败，忽略或抛出异常
                    System.err.println("无法解析 SRA 编号范围: " + part);
                }
            } else {
                // 单个 SRA 编号
                result.add(part);
            }
        }
        return result.toArray(new String[0]);
    }

    private static void writeString2Txt(String context) {
        try (FileWriter fileWriter = new FileWriter("./sra.txt")) {
            fileWriter.write(context);
        } catch (IOException e) {
            throw new RuntimeException("写入 sra.txt 文件失败", e);
        }
    }
}