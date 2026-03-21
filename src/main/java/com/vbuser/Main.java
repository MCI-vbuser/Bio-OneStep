package com.vbuser;

import com.vbuser.fastq.FastqProcessor;
import com.vbuser.fastq.SRA2FastQ;
import com.vbuser.pre.CommonEnv;
import com.vbuser.pre.DownloadSRA;

import java.io.File;
import java.io.IOException;

public class Main {
    public static void main(String[] args){
        try {
            WebConsole.start();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        CommonEnv.envCheck();
        cls();
        System.out.println("环境检查通过");
        prepareSRA();
        WebConsole.stopServer();
        try {
            Thread.sleep(1000);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
        WebConsole.killProcess();
    }

    private static void cls(){
        try {
            new ProcessBuilder("bash", "-c", "clear").inheritIO().start().waitFor();
        } catch (InterruptedException | IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void prepareSRA(){
        System.out.println("[1]使用本地SRA文件\n[2]下载新的SRA文件\n[3]已有bam,跳过");
        String ans = WebConsole.readLine();
        if(ans.equals("1")){
            System.out.println("请指定SRA文件路径");
        }else if(ans.equals("2")){
            System.out.println("请指定下载的SRR编号");
            String srr = WebConsole.readLine();
            DownloadSRA.download(srr);
        } else{
            cls();
            return;
        }
        prepareFastQ(Integer.parseInt(ans));
    }

    private static void prepareFastQ(int choice){
        if (choice == 2) {
            SRA2FastQ.handle(new File("./sra"));
        } else if (choice == 1) {
            SRA2FastQ.handle(new File(WebConsole.readLine()));
        } else return;
        System.out.println("fastq文件已就位，开始组装bam文件");
        prepareBam();
        cls();
        System.out.println("组装bam完成");
    }

    private static void prepareBam(){
        try {
            FastqProcessor.processFastqToBam();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

}
