package com.vbuser;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import com.sun.management.OperatingSystemMXBean;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileStore;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class WebConsole {

    private static final int PORT = 11451;
    private static HttpServer server;  // 保存服务器引用，用于关闭

    public static void start() throws IOException {
        // 1. 重定向 System.out 到网页 + 原始控制台
        OutputCapture.init();
        // 2. 重定向 System.in 并初始化输入读取器
        InputRedirect.init();
        // 3. 启动系统监控后台线程
        SystemMonitor.start();

        // 4. 启动 HTTP 服务器（设置守护线程池）
        server = HttpServer.create(new InetSocketAddress(PORT), 0);
        server.createContext("/", new WebHandler());
        server.createContext("/logs", new LogsHandler());
        server.createContext("/input", new InputHandler());
        server.createContext("/stats", new StatsHandler());
        server.createContext("/stop", new StopHandler());  // 新增终止端点

        // 使用守护线程池，主线程结束后 JVM 自动退出
        server.setExecutor(Executors.newCachedThreadPool(r -> {
            Thread t = new Thread(r);
            t.setDaemon(true);
            return t;
        }));
        server.start();
        System.out.println("Web Console started at http://localhost:" + PORT);
    }

    /**
     * 终止整个进程（强制退出）
     */
    public static void killProcess() {
        System.exit(0);
    }

    /**
     * 停止 HTTP 服务器（不退出 JVM）
     */
    public static void stopServer() {
        if (server != null) {
            server.stop(0);
        }
    }

    /**
     * 替代 System.console().readLine() 的方法
     * @param prompt 提示信息（会输出到网页）
     * @return 用户输入的一行字符串（不含换行符）
     */
    public static String readLine(String prompt) {
        System.out.print(prompt);
        try {
            return inputReader.readLine();
        } catch (IOException e) {
            return null;
        }
    }

    public static String readLine() {
        return readLine("");
    }

    private static BufferedReader inputReader; // 供 readLine 使用

    // ========== 输出捕获 ==========
    static class OutputCapture extends PrintStream {
        private static final List<String> logList = new ArrayList<>();
        private static final ReentrantReadWriteLock lock = new ReentrantReadWriteLock();
        private static final PrintStream originalOut = System.out;

        private OutputCapture(OutputStream out) {
            super(out);
        }

        static void init() {
            System.setOut(new OutputCapture(new FilterOutputStream(originalOut) {
                @Override
                public void write(byte[] b, int off, int len) {
                    originalOut.write(b, off, len);
                    String line = new String(b, off, len, StandardCharsets.UTF_8);
                    lock.writeLock().lock();
                    try {
                        logList.add(line);
                    } finally {
                        lock.writeLock().unlock();
                    }
                }

                @Override
                public void write(int b) {
                    originalOut.write(b);
                }
            }));
        }

        static LogResponse getLogsSince(int lastIndex) {
            lock.readLock().lock();
            try {
                int size = logList.size();
                if (lastIndex >= size - 1) {
                    return new LogResponse(null, size - 1);
                }
                int start = Math.max(0, lastIndex + 1);
                List<String> newLogs = new ArrayList<>(logList.subList(start, size));
                return new LogResponse(newLogs, size - 1);
            } finally {
                lock.readLock().unlock();
            }
        }
    }

    static class LogResponse {
        final List<String> logs;
        final int lastIndex;
        LogResponse(List<String> logs, int lastIndex) {
            this.logs = logs;
            this.lastIndex = lastIndex;
        }
    }

    // ========== 输入重定向 ==========
    static class InputRedirect {
        // 队列容量根据实际需求调整，此处设为 1000
        private static final int QUEUE_CAPACITY = 1000;
        private static final BlockingQueue<byte[]> inputQueue = new LinkedBlockingQueue<>(QUEUE_CAPACITY);

        /**
         * 初始化输入重定向，将 System.in 替换为自定义的 QueueInputStream
         */
        static void init() {
            System.setIn(new QueueInputStream(inputQueue));
            // 供 WebConsole.readLine 使用的 BufferedReader
            WebConsole.inputReader = new BufferedReader(new InputStreamReader(System.in, StandardCharsets.UTF_8));
        }

        /**
         * 发送用户输入（由 HTTP 请求调用）
         * @param line 用户输入的一行字符串（不含换行符）
         * @throws IOException 当队列满时抛出，HTTP 层可据此返回 500 错误
         */
        static void sendInput(String line) throws IOException {
            byte[] data = (line + "\n").getBytes(StandardCharsets.UTF_8);
            boolean offered = inputQueue.offer(data);
            if (!offered) {
                throw new IOException("Input queue is full, please try again later.");
            }
        }

        /**
         * 自定义 InputStream，从 BlockingQueue 中读取字节数组
         */
        static class QueueInputStream extends InputStream {
            private final BlockingQueue<byte[]> queue;
            private byte[] currentBuffer;
            private int pos;

            QueueInputStream(BlockingQueue<byte[]> queue) {
                this.queue = queue;
            }

            @Override
            public int read() throws IOException {
                if (currentBuffer == null || pos >= currentBuffer.length) {
                    try {
                        currentBuffer = queue.take();  // 阻塞等待输入
                        pos = 0;
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        throw new IOException(e);
                    }
                }
                return currentBuffer[pos++] & 0xFF;
            }

            @Override
            public int read(byte[] b, int off, int len) throws IOException {
                if (len == 0) return 0;
                int c = read();
                if (c == -1) return -1;
                b[off] = (byte) c;
                return 1;
            }
        }
    }

    // ========== 系统监控 ==========
    static class SystemMonitor {
        private static volatile StatsData currentStats = new StatsData(0, 0, 0, 0, 0, 0, 0, 0, 0);
        private static final Object lock = new Object();
        private static final String networkInterface; // 可配置的网络接口

        static {
            // 从系统属性读取网络接口名，默认 null 表示自动选择
            networkInterface = System.getProperty("webconsole.network.interface");
        }

        static void start() {
            Thread monitorThread = new Thread(() -> {
                long lastNetTime = System.currentTimeMillis();
                long lastRx = 0, lastTx = 0;
                // 确定网络接口
                String iface = determineNetworkInterface();
                if (iface == null) {
                    System.err.println("[WebConsole] 无法确定网络接口，网络速率将显示为0");
                } else {
                    System.out.println("[WebConsole] 监控网络接口: " + iface);
                }

                while (true) {
                    try {
                        Thread.sleep(1000); // 每秒采集一次
                        long now = System.currentTimeMillis();
                        double elapsedSec = (now - lastNetTime) / 1000.0;

                        // CPU 使用率
                        double cpuUsage = getCpuUsage();

                        // 内存信息（改进版：模拟 MemAvailable）
                        MemInfo mem = getMemInfo();
                        long totalMem = mem.total;
                        long usedMem = mem.used;
                        double memUsagePercent = totalMem > 0 ? (usedMem * 100.0 / totalMem) : 0;

                        // 磁盘信息（根目录）
                        DiskInfo disk = getDiskInfo();

                        // 网络信息（速率）
                        long rxBytes = 0, txBytes = 0;
                        if (iface != null) {
                            NetworkStats net = getNetworkStats(iface);
                            rxBytes = net.rxBytes;
                            txBytes = net.txBytes;
                        }
                        double rxRate = (rxBytes - lastRx) / elapsedSec;
                        double txRate = (txBytes - lastTx) / elapsedSec;
                        if (rxRate < 0) rxRate = 0;
                        if (txRate < 0) txRate = 0;
                        lastRx = rxBytes;
                        lastTx = txBytes;
                        lastNetTime = now;

                        // 更新当前数据
                        synchronized (lock) {
                            currentStats = new StatsData(
                                    cpuUsage,
                                    totalMem, usedMem, memUsagePercent,
                                    disk.total, disk.used, disk.usagePercent,
                                    rxRate, txRate
                            );
                        }
                    } catch (InterruptedException e) {
                        break; // 线程被中断，退出
                    } catch (Exception e) {
                        // 其他异常不中断线程，打印错误后继续
                        throw new RuntimeException(e);
                    }
                }
            });
            monitorThread.setDaemon(true);
            monitorThread.start();
        }

        static StatsData getStats() {
            synchronized (lock) {
                return currentStats;
            }
        }

        private static double getCpuUsage() {
            try {
                OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
                double load = osBean.getSystemCpuLoad();
                return load >= 0 ? load * 100 : -1;
            } catch (Exception e) {
                return -1;
            }
        }

        /**
         * 获取内存信息（总内存、已用内存）
         * 优先从 /proc/meminfo 获取 MemTotal 和 (MemAvailable 或 MemFree+Buffers+Cached)
         */
        private static MemInfo getMemInfo() {
            long total = 0;
            long available = 0;
            try (BufferedReader reader = new BufferedReader(new FileReader("/proc/meminfo"))) {
                String line;
                long memFree = 0, buffers = 0, cached = 0;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("MemTotal:")) {
                        total = parseMeminfoValue(line);
                    } else if (line.startsWith("MemAvailable:")) {
                        available = parseMeminfoValue(line);
                    } else if (line.startsWith("MemFree:")) {
                        memFree = parseMeminfoValue(line);
                    } else if (line.startsWith("Buffers:")) {
                        buffers = parseMeminfoValue(line);
                    } else if (line.startsWith("Cached:")) {
                        cached = parseMeminfoValue(line);
                    }
                }
                // 优先使用 MemAvailable（内核2.6.27+）
                if (total > 0 && available > 0) {
                    long used = total - available;
                    return new MemInfo(total, used);
                }
                // 降级：模拟 MemAvailable = MemFree + Buffers + Cached
                if (total > 0) {
                    long simulatedAvailable = memFree + buffers + cached;
                    if (simulatedAvailable > 0 && simulatedAvailable <= total) {
                        long used = total - simulatedAvailable;
                        return new MemInfo(total, used);
                    }
                }
            } catch (Exception e) {
                // 读取失败，使用 MXBean
            }

            // 最终回退到 MXBean
            try {
                OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
                total = osBean.getTotalPhysicalMemorySize();
                long free = osBean.getFreePhysicalMemorySize();
                long used = total - free;
                return new MemInfo(total, used);
            } catch (Exception e) {
                return new MemInfo(0, 0);
            }
        }

        private static long parseMeminfoValue(String line) {
            // 格式: "MemTotal:       16384000 kB"
            String[] parts = line.trim().split("\\s+");
            if (parts.length >= 2) {
                try {
                    long value = Long.parseLong(parts[1]);
                    // 如果单位是 kB，转换为字节
                    if (parts.length >= 3 && "kB".equalsIgnoreCase(parts[2])) {
                        value *= 1024;
                    }
                    return value;
                } catch (NumberFormatException e) {
                    return 0;
                }
            }
            return 0;
        }

        private static DiskInfo getDiskInfo() {
            try {
                FileStore store = Files.getFileStore(Paths.get("/"));
                long total = store.getTotalSpace();
                long free = store.getUsableSpace();
                long used = total - free;
                double percent = total > 0 ? (used * 100.0 / total) : 0;
                return new DiskInfo(total, used, percent);
            } catch (Exception e) {
                return new DiskInfo(0, 0, 0);
            }
        }

        /**
         * 确定要监控的网络接口
         * 1. 若指定了系统属性 webconsole.network.interface，使用该接口
         * 2. 否则从 /proc/net/dev 中找出第一个非 lo 且接收/发送字节数>0的接口
         */
        private static String determineNetworkInterface() {
            if (networkInterface != null && !networkInterface.isEmpty()) {
                return networkInterface;
            }

            try (BufferedReader reader = new BufferedReader(new FileReader("/proc/net/dev"))) {
                String line;
                boolean firstLine = true;
                while ((line = reader.readLine()) != null) {
                    if (firstLine) {
                        firstLine = false;
                        continue;
                    }
                    line = line.trim();
                    if (line.isEmpty()) continue;
                    String[] parts = line.split("\\s+");
                    if (parts.length < 10) continue;
                    String iface = parts[0].replace(":", "");
                    if (iface.equals("lo")) continue;
                    // 尝试解析接收和发送字节数，如果不为0则优先选择有流量的接口
                    try {
                        long rx = Long.parseLong(parts[1]);
                        long tx = Long.parseLong(parts[9]);
                        if (rx > 0 || tx > 0) {
                            return iface;
                        }
                    } catch (NumberFormatException e) {
                        // 忽略，继续下一个接口
                    }
                }
                // 如果没有找到有流量的接口，返回第一个非 lo 的接口
                try (BufferedReader reader2 = new BufferedReader(new FileReader("/proc/net/dev"))) {
                    boolean first = true;
                    while ((line = reader2.readLine()) != null) {
                        if (first) {
                            first = false;
                            continue;
                        }
                        line = line.trim();
                        if (line.isEmpty()) continue;
                        String[] parts = line.split("\\s+");
                        if (parts.length >= 1) {
                            String iface = parts[0].replace(":", "");
                            if (!iface.equals("lo")) {
                                return iface;
                            }
                        }
                    }
                }
            } catch (Exception e) {
                // 忽略
            }
            return null;
        }

        private static NetworkStats getNetworkStats(String iface) {
            long rx = 0, tx = 0;
            try (BufferedReader reader = new BufferedReader(new FileReader("/proc/net/dev"))) {
                String line;
                boolean firstLine = true;
                while ((line = reader.readLine()) != null) {
                    if (firstLine) {
                        firstLine = false;
                        continue;
                    }
                    line = line.trim();
                    if (line.isEmpty()) continue;
                    String[] parts = line.split("\\s+");
                    if (parts.length < 10) continue;
                    String name = parts[0].replace(":", "");
                    if (name.equals(iface)) {
                        rx = Long.parseLong(parts[1]);
                        tx = Long.parseLong(parts[9]);
                        break;
                    }
                }
            } catch (Exception e) {
                // 忽略
            }
            return new NetworkStats(rx, tx);
        }

        static class StatsData {
            final double cpuUsage;
            final long totalMem;
            final long usedMem;
            final double memUsagePercent;
            final long totalDisk;
            final long usedDisk;
            final double diskUsagePercent;
            final double rxRate;  // bytes/sec
            final double txRate;

            StatsData(double cpuUsage, long totalMem, long usedMem, double memUsagePercent,
                      long totalDisk, long usedDisk, double diskUsagePercent,
                      double rxRate, double txRate) {
                this.cpuUsage = cpuUsage;
                this.totalMem = totalMem;
                this.usedMem = usedMem;
                this.memUsagePercent = memUsagePercent;
                this.totalDisk = totalDisk;
                this.usedDisk = usedDisk;
                this.diskUsagePercent = diskUsagePercent;
                this.rxRate = rxRate;
                this.txRate = txRate;
            }
        }

        static class MemInfo {
            long total, used;
            MemInfo(long total, long used) {
                this.total = total;
                this.used = used;
            }
        }

        static class DiskInfo {
            long total, used;
            double usagePercent;
            DiskInfo(long total, long used, double usagePercent) {
                this.total = total;
                this.used = used;
                this.usagePercent = usagePercent;
            }
        }

        static class NetworkStats {
            long rxBytes, txBytes;
            NetworkStats(long rx, long tx) {
                this.rxBytes = rx;
                this.txBytes = tx;
            }
        }
    }

    // ========== HTTP 处理器 ==========
    static class WebHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            String html = getHtmlPage();
            exchange.getResponseHeaders().set("Content-Type", "text/html; charset=UTF-8");
            exchange.sendResponseHeaders(200, html.getBytes(StandardCharsets.UTF_8).length);
            try (OutputStream os = exchange.getResponseBody()) {
                os.write(html.getBytes(StandardCharsets.UTF_8));
            }
        }

        private String getHtmlPage() {
            return "<!DOCTYPE html>\n" +
                    "<html lang=\"zh-CN\">\n" +
                    "<head>\n" +
                    "    <meta charset=\"UTF-8\">\n" +
                    "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n" +
                    "    <title>网页控制台</title>\n" +
                    "    <style>\n" +
                    "        body {\n" +
                    "            background-color: #0d1117;\n" +
                    "            color: #c9d1d9;\n" +
                    "            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Noto Sans', Helvetica, Arial, sans-serif;\n" +
                    "            margin: 0;\n" +
                    "            padding: 20px;\n" +
                    "        }\n" +
                    "        .container {\n" +
                    "            max-width: 1400px;\n" +
                    "            margin: 0 auto;\n" +
                    "        }\n" +
                    "        .header {\n" +
                    "            border-bottom: 1px solid #30363d;\n" +
                    "            padding-bottom: 10px;\n" +
                    "            margin-bottom: 20px;\n" +
                    "        }\n" +
                    "        .header h1 {\n" +
                    "            font-size: 24px;\n" +
                    "            font-weight: 600;\n" +
                    "            margin: 0;\n" +
                    "        }\n" +
                    "        /* 监控卡片区域 */\n" +
                    "        .stats-grid {\n" +
                    "            display: grid;\n" +
                    "            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));\n" +
                    "            gap: 16px;\n" +
                    "            margin-bottom: 24px;\n" +
                    "        }\n" +
                    "        .stat-card {\n" +
                    "            background-color: #161b22;\n" +
                    "            border: 1px solid #30363d;\n" +
                    "            border-radius: 6px;\n" +
                    "            padding: 12px 16px;\n" +
                    "        }\n" +
                    "        .stat-title {\n" +
                    "            font-size: 12px;\n" +
                    "            font-weight: 500;\n" +
                    "            color: #8b949e;\n" +
                    "            margin-bottom: 8px;\n" +
                    "            text-transform: uppercase;\n" +
                    "            letter-spacing: 0.5px;\n" +
                    "        }\n" +
                    "        .stat-value {\n" +
                    "            font-size: 24px;\n" +
                    "            font-weight: 600;\n" +
                    "            margin-bottom: 8px;\n" +
                    "        }\n" +
                    "        .stat-sub {\n" +
                    "            font-size: 12px;\n" +
                    "            color: #8b949e;\n" +
                    "        }\n" +
                    "        .progress-bar {\n" +
                    "            background-color: #30363d;\n" +
                    "            border-radius: 10px;\n" +
                    "            height: 6px;\n" +
                    "            overflow: hidden;\n" +
                    "            margin-top: 8px;\n" +
                    "        }\n" +
                    "        .progress-fill {\n" +
                    "            background-color: #238636;\n" +
                    "            height: 100%;\n" +
                    "            width: 0%;\n" +
                    "            border-radius: 10px;\n" +
                    "        }\n" +
                    "        /* 控制台区域 */\n" +
                    "        .console {\n" +
                    "            background-color: #161b22;\n" +
                    "            border: 1px solid #30363d;\n" +
                    "            border-radius: 6px;\n" +
                    "            padding: 16px;\n" +
                    "            font-family: 'SFMono', 'Monaco', 'Cascadia Code', monospace;\n" +
                    "            font-size: 12px;\n" +
                    "            line-height: 1.5;\n" +
                    "            overflow-y: auto;\n" +
                    "            height: 400px;\n" +
                    "            margin-bottom: 20px;\n" +
                    "        }\n" +
                    "        .console-line {\n" +
                    "            white-space: pre-wrap;\n" +
                    "            word-break: break-all;\n" +
                    "            border-bottom: 1px solid #21262d;\n" +
                    "            padding: 2px 0;\n" +
                    "        }\n" +
                    "        .input-area {\n" +
                    "            display: flex;\n" +
                    "            gap: 10px;\n" +
                    "            align-items: center;\n" +
                    "        }\n" +
                    "        .input-area input {\n" +
                    "            flex: 1;\n" +
                    "            background-color: #0d1117;\n" +
                    "            border: 1px solid #30363d;\n" +
                    "            border-radius: 6px;\n" +
                    "            color: #c9d1d9;\n" +
                    "            padding: 8px 12px;\n" +
                    "            font-family: inherit;\n" +
                    "            font-size: 14px;\n" +
                    "        }\n" +
                    "        .input-area button {\n" +
                    "            background-color: #238636;\n" +
                    "            border: none;\n" +
                    "            border-radius: 6px;\n" +
                    "            color: white;\n" +
                    "            padding: 8px 16px;\n" +
                    "            font-size: 14px;\n" +
                    "            font-weight: 500;\n" +
                    "            cursor: pointer;\n" +
                    "            transition: background-color 0.2s;\n" +
                    "        }\n" +
                    "        .input-area button:hover {\n" +
                    "            background-color: #2ea043;\n" +
                    "        }\n" +
                    "        .kill-btn {\n" +
                    "            background-color: #da3633 !important;\n" +
                    "        }\n" +
                    "        .kill-btn:hover {\n" +
                    "            background-color: #f85149 !important;\n" +
                    "        }\n" +
                    "        .status {\n" +
                    "            margin-top: 10px;\n" +
                    "            font-size: 12px;\n" +
                    "            color: #8b949e;\n" +
                    "        }\n" +
                    "    </style>\n" +
                    "</head>\n" +
                    "<body>\n" +
                    "<div class=\"container\">\n" +
                    "    <div class=\"header\">\n" +
                    "        <h1>实时控制台输出 + 系统监控</h1>\n" +
                    "    </div>\n" +
                    "\n" +
                    "    <!-- 监控卡片区域 -->\n" +
                    "    <div class=\"stats-grid\" id=\"statsGrid\">\n" +
                    "        <div class=\"stat-card\">\n" +
                    "            <div class=\"stat-title\">CPU 占用</div>\n" +
                    "            <div class=\"stat-value\" id=\"cpuValue\">--%</div>\n" +
                    "            <div class=\"progress-bar\"><div class=\"progress-fill\" id=\"cpuFill\"></div></div>\n" +
                    "        </div>\n" +
                    "        <div class=\"stat-card\">\n" +
                    "            <div class=\"stat-title\">内存占用</div>\n" +
                    "            <div class=\"stat-value\" id=\"memValue\">-- / --</div>\n" +
                    "            <div class=\"progress-bar\"><div class=\"progress-fill\" id=\"memFill\"></div></div>\n" +
                    "            <div class=\"stat-sub\" id=\"memPercent\">--%</div>\n" +
                    "        </div>\n" +
                    "        <div class=\"stat-card\">\n" +
                    "            <div class=\"stat-title\">磁盘占用</div>\n" +
                    "            <div class=\"stat-value\" id=\"diskValue\">-- / --</div>\n" +
                    "            <div class=\"progress-bar\"><div class=\"progress-fill\" id=\"diskFill\"></div></div>\n" +
                    "            <div class=\"stat-sub\" id=\"diskPercent\">--%</div>\n" +
                    "        </div>\n" +
                    "        <div class=\"stat-card\">\n" +
                    "            <div class=\"stat-title\">网络速率</div>\n" +
                    "            <div class=\"stat-value\" id=\"netRx\">↓ -- KB/s</div>\n" +
                    "            <div class=\"stat-value\" id=\"netTx\" style=\"margin-top: 4px;\">↑ -- KB/s</div>\n" +
                    "        </div>\n" +
                    "    </div>\n" +
                    "\n" +
                    "    <!-- 控制台区域 -->\n" +
                    "    <div class=\"console\" id=\"console\">\n" +
                    "        <div class=\"console-line\">等待日志输出...</div>\n" +
                    "    </div>\n" +
                    "    <div class=\"input-area\">\n" +
                    "        <input type=\"text\" id=\"inputField\" placeholder=\"输入命令并回车或点击发送...\" autocomplete=\"off\">\n" +
                    "        <button id=\"sendBtn\">发送</button>\n" +
                    "        <button id=\"killBtn\" class=\"kill-btn\">自刎归天</button>\n" +
                    "    </div>\n" +
                    "    <div class=\"status\" id=\"status\">已连接</div>\n" +
                    "</div>\n" +
                    "\n" +
                    "<script>\n" +
                    "    // 日志轮询\n" +
                    "    let lastIndex = -1;\n" +
                    "    const consoleDiv = document.getElementById('console');\n" +
                    "    const inputField = document.getElementById('inputField');\n" +
                    "    const sendBtn = document.getElementById('sendBtn');\n" +
                    "    const statusSpan = document.getElementById('status');\n" +
                    "    const killBtn = document.getElementById('killBtn');\n" +
                    "\n" +
                    "    function fetchLogs() {\n" +
                    "        fetch('/logs?last=' + lastIndex)\n" +
                    "            .then(response => response.json())\n" +
                    "            .then(data => {\n" +
                    "                if (data.logs && data.logs.length > 0) {\n" +
                    "                    data.logs.forEach(line => {\n" +
                    "                        const lineDiv = document.createElement('div');\n" +
                    "                        lineDiv.className = 'console-line';\n" +
                    "                        lineDiv.textContent = line;\n" +
                    "                        consoleDiv.appendChild(lineDiv);\n" +
                    "                    });\n" +
                    "                    consoleDiv.scrollTop = consoleDiv.scrollHeight;\n" +
                    "                    lastIndex = data.lastIndex;\n" +
                    "                }\n" +
                    "                statusSpan.textContent = '已连接';\n" +
                    "            })\n" +
                    "            .catch(err => {\n" +
                    "                statusSpan.textContent = '连接错误: ' + err.message;\n" +
                    "            });\n" +
                    "    }\n" +
                    "\n" +
                    "    function sendInput() {\n" +
                    "        const input = inputField.value;\n" +
                    "        if (!input.trim()) return;\n" +
                    "        fetch('/input', {\n" +
                    "            method: 'POST',\n" +
                    "            body: input\n" +
                    "        })\n" +
                    "        .then(response => {\n" +
                    "            if (response.ok) {\n" +
                    "                inputField.value = '';\n" +
                    "                statusSpan.textContent = '输入已发送';\n" +
                    "            } else {\n" +
                    "                statusSpan.textContent = '发送失败';\n" +
                    "            }\n" +
                    "        })\n" +
                    "        .catch(err => {\n" +
                    "            statusSpan.textContent = '发送错误: ' + err.message;\n" +
                    "        });\n" +
                    "    }\n" +
                    "\n" +
                    "    // 监控数据轮询\n" +
                    "    function fetchStats() {\n" +
                    "        fetch('/stats')\n" +
                    "            .then(response => response.json())\n" +
                    "            .then(data => {\n" +
                    "                // CPU\n" +
                    "                const cpu = data.cpuUsage;\n" +
                    "                const cpuVal = cpu >= 0 ? cpu.toFixed(1) + '%' : 'N/A';\n" +
                    "                document.getElementById('cpuValue').innerText = cpuVal;\n" +
                    "                document.getElementById('cpuFill').style.width = (cpu >=0 ? cpu : 0) + '%';\n" +
                    "\n" +
                    "                // 内存\n" +
                    "                const memTotal = (data.totalMem / (1024**3)).toFixed(1);\n" +
                    "                const memUsed = (data.usedMem / (1024**3)).toFixed(1);\n" +
                    "                document.getElementById('memValue').innerText = memUsed + ' GB / ' + memTotal + ' GB';\n" +
                    "                const memPercent = data.memUsagePercent.toFixed(1);\n" +
                    "                document.getElementById('memPercent').innerText = memPercent + '%';\n" +
                    "                document.getElementById('memFill').style.width = memPercent + '%';\n" +
                    "\n" +
                    "                // 磁盘\n" +
                    "                const diskTotal = (data.totalDisk / (1024**3)).toFixed(1);\n" +
                    "                const diskUsed = (data.usedDisk / (1024**3)).toFixed(1);\n" +
                    "                document.getElementById('diskValue').innerText = diskUsed + ' GB / ' + diskTotal + ' GB';\n" +
                    "                const diskPercent = data.diskUsagePercent.toFixed(1);\n" +
                    "                document.getElementById('diskPercent').innerText = diskPercent + '%';\n" +
                    "                document.getElementById('diskFill').style.width = diskPercent + '%';\n" +
                    "\n" +
                    "                // 网络\n" +
                    "                const rxKbps = data.rxRate / 1024;\n" +
                    "                const txKbps = data.txRate / 1024;\n" +
                    "                document.getElementById('netRx').innerHTML = '↓ ' + rxKbps.toFixed(1) + ' KB/s';\n" +
                    "                document.getElementById('netTx').innerHTML = '↑ ' + txKbps.toFixed(1) + ' KB/s';\n" +
                    "            })\n" +
                    "            .catch(err => console.error('Stats fetch error:', err));\n" +
                    "    }\n" +
                    "\n" +
                    "    // 终止进程\n" +
                    "    killBtn.addEventListener('click', () => {\n" +
                    "        if (confirm('确定要终止整个进程吗？此操作不可逆。')) {\n" +
                    "            fetch('/stop', { method: 'POST' }).then(() => {\n" +
                    "                statusSpan.textContent = '进程终止中...';\n" +
                    "            }).catch(() => {\n" +
                    "                statusSpan.textContent = '进程已终止';\n" +
                    "            });\n" +
                    "        }\n" +
                    "    });\n" +
                    "\n" +
                    "    inputField.addEventListener('keypress', (e) => {\n" +
                    "        if (e.key === 'Enter') sendInput();\n" +
                    "    });\n" +
                    "    sendBtn.addEventListener('click', sendInput);\n" +
                    "\n" +
                    "    setInterval(fetchLogs, 500);\n" +
                    "    fetchLogs();\n" +
                    "    setInterval(fetchStats, 2000);\n" +
                    "    fetchStats();\n" +
                    "</script>\n" +
                    "</body>\n" +
                    "</html>";
        }
    }

    static class LogsHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            String query = exchange.getRequestURI().getQuery();
            int lastIndex = -1;
            if (query != null && query.startsWith("last=")) {
                try {
                    lastIndex = Integer.parseInt(query.substring(5));
                } catch (NumberFormatException ignored) {}
            }
            LogResponse resp = OutputCapture.getLogsSince(lastIndex);
            String json;
            if (resp.logs == null || resp.logs.isEmpty()) {
                json = "{\"logs\":[],\"lastIndex\":" + resp.lastIndex + "}";
            } else {
                StringBuilder sb = new StringBuilder();
                sb.append("{\"logs\":[");
                for (int i = 0; i < resp.logs.size(); i++) {
                    if (i > 0) sb.append(",");
                    sb.append("\"").append(escapeJson(resp.logs.get(i))).append("\"");
                }
                sb.append("],\"lastIndex\":").append(resp.lastIndex).append("}");
                json = sb.toString();
            }
            exchange.getResponseHeaders().set("Content-Type", "application/json; charset=UTF-8");
            exchange.sendResponseHeaders(200, json.getBytes(StandardCharsets.UTF_8).length);
            try (OutputStream os = exchange.getResponseBody()) {
                os.write(json.getBytes(StandardCharsets.UTF_8));
            }
        }

        private String escapeJson(String s) {
            return s.replace("\\", "\\\\")
                    .replace("\"", "\\\"")
                    .replace("\n", "\\n")
                    .replace("\r", "\\r");
        }
    }

    static class InputHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            if (!"POST".equalsIgnoreCase(exchange.getRequestMethod())) {
                exchange.sendResponseHeaders(405, -1);
                return;
            }
            String body;
            try (InputStream is = exchange.getRequestBody();
                 BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8))) {
                StringBuilder sb = new StringBuilder();
                String line;
                while ((line = reader.readLine()) != null) {
                    sb.append(line);
                }
                body = sb.toString();
            }
            try {
                InputRedirect.sendInput(body);
                exchange.sendResponseHeaders(200, 0);
            } catch (Exception e) {
                exchange.sendResponseHeaders(500, 0);
            } finally {
                exchange.close();
            }
        }
    }

    static class StatsHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            String json = getString();
            exchange.getResponseHeaders().set("Content-Type", "application/json; charset=UTF-8");
            exchange.sendResponseHeaders(200, json.getBytes(StandardCharsets.UTF_8).length);
            try (OutputStream os = exchange.getResponseBody()) {
                os.write(json.getBytes(StandardCharsets.UTF_8));
            }
        }

        private String getString() {
            SystemMonitor.StatsData stats = SystemMonitor.getStats();
            return String.format(Locale.US,
                    "{\"cpuUsage\":%.2f,\"totalMem\":%d,\"usedMem\":%d,\"memUsagePercent\":%.2f,\"totalDisk\":%d,\"usedDisk\":%d,\"diskUsagePercent\":%.2f,\"rxRate\":%.2f,\"txRate\":%.2f}",
                    stats.cpuUsage,
                    stats.totalMem, stats.usedMem, stats.memUsagePercent,
                    stats.totalDisk, stats.usedDisk, stats.diskUsagePercent,
                    stats.rxRate, stats.txRate);
        }
    }

    // ========== 新增 /stop 处理器 ==========
    static class StopHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            if (!"POST".equalsIgnoreCase(exchange.getRequestMethod())) {
                exchange.sendResponseHeaders(405, -1);
                return;
            }

            String response = "Process is shutting down...";
            exchange.getResponseHeaders().set("Content-Type", "text/plain; charset=UTF-8");
            exchange.sendResponseHeaders(200, response.getBytes(StandardCharsets.UTF_8).length);
            try (OutputStream os = exchange.getResponseBody()) {
                os.write(response.getBytes(StandardCharsets.UTF_8));
            }

            // 延迟退出，确保响应已发送
            new Thread(() -> {
                try {
                    Thread.sleep(100);
                } catch (InterruptedException ignored) {}
                System.exit(0);
            }).start();
        }
    }
}