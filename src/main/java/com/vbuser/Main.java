package com.vbuser;

import com.vbuser.pre.CommonEnv;
import com.vbuser.pre.DownloadSRA;

import java.io.IOException;

public class Main {
    public static void main(String[] args){
        CommonEnv.envCheck();
        cls();
        System.out.println("环境检查通过");
        prepareSRA();
    }

    private static void cls(){
        try {
            new ProcessBuilder("bash", "-c", "clear").inheritIO().start().waitFor();
        } catch (InterruptedException | IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void prepareSRA(){
        System.out.println("[1]使用本地SRA文件\n[2]下载新的SRA文件\n[3]跳过");
        String ans = System.console().readLine();
        if(ans.equals("1")){
            System.out.println("请指定SRA文件路径");
        }else if(ans.equals("2")){
            System.out.println("请指定下载的SRR编号");
            String srr = System.console().readLine();
            DownloadSRA.download(srr);
        } else{
            cls();
        }
    }

}
