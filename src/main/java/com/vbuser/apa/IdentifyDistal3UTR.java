package com.vbuser.apa;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * High-performance Java reimplementation of identifyDistal3UTR.pl.
 * Strictly algorithm-equivalent to Perl version, optimized with prefix sums and parallel processing.
 * Updated: double precision, edge-case alignment for identical output.
 */
public class IdentifyDistal3UTR {

    // Parameters
    private static int windowSize = 100;
    private static int extendSize = 10000;
    private static double coverageCutoff = 0.05;
    private static double percentageCutoff = 0.8;
    private static String geneSymbolFile = null;
    private static List<String> inputFiles;
    private static String geneModelFile;
    private static String outputFile;
    private static int numThreads = Runtime.getRuntime().availableProcessors();

    // Data structures
    private static final Map<String, UtrInfo> utrMap = new LinkedHashMap<>();
    private static final Map<String, Map<String, UtrInfo>> chrUtrMap = new HashMap<>();
    private static final List<SampleData> samples = new ArrayList<>();

    // Statistics & progress
    private static final AtomicInteger processedUtrs = new AtomicInteger(0);
    private static volatile boolean done = false;
    private static long startTime;

    public static void main(String[] args) throws Exception {
        parseArguments(args);
        validateArguments();

        System.err.println("Extracting annotated 3'UTRs...");
        loadGeneModel();

        System.err.println("Loading bedgraph files...");
        loadBedgraphFiles();

        System.err.println("Identifying distal 3'UTRs...");
        identifyDistalUTRs();

        System.err.println("Removing overlapping UTRs...");
        removeOverlaps();

        System.err.println("Writing output...");
        writeOutput();

        System.err.println("Done.");
    }

    private static void parseArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-i":
                    inputFiles = new ArrayList<>();
                    while (i + 1 < args.length && !args[i + 1].startsWith("-")) inputFiles.add(args[++i]);
                    break;
                case "-m": geneModelFile = args[++i]; break;
                case "-o": outputFile = args[++i]; break;
                case "-w": windowSize = Integer.parseInt(args[++i]); break;
                case "-e": extendSize = Integer.parseInt(args[++i]); break;
                case "-c": coverageCutoff = Double.parseDouble(args[++i]); break;
                case "-p": percentageCutoff = Double.parseDouble(args[++i]); break;
                case "-s": geneSymbolFile = args[++i]; break;
                case "-t": numThreads = Integer.parseInt(args[++i]); break;
                case "-h": printUsage(); System.exit(0);
            }
        }
    }

    private static void printUsage() {
        System.err.println("Usage: java IdentifyDistal3UTR -i <bedgraph...> -m <genemodel.bed> -o <output.bed> [options]");
        System.err.println("Options: -w <int> -e <int> -c <double> -p <double> -s <symbols.txt> -t <threads>");
    }

    private static void validateArguments() {
        if (inputFiles == null || inputFiles.isEmpty()) error("No input files (-i)");
        if (geneModelFile == null) error("Gene model file (-m) required");
        if (outputFile == null) error("Output file (-o) required");
        if (windowSize < 20) error("Window size (-w) must be >= 20");
        if (extendSize < 100) error("Extend size (-e) must be >= 100");
        if (coverageCutoff <= 0 || coverageCutoff >= 1) error("Coverage cutoff (-c) in (0,1)");
        if (percentageCutoff <= 0 || percentageCutoff >= 1) error("Percentage cutoff (-p) in (0,1)");
        if (geneSymbolFile != null && !Files.exists(Paths.get(geneSymbolFile))) error("Symbol file not found");
    }

    private static void error(String msg) {
        throw new IllegalArgumentException("Error: " + msg);
    }

    // --------------------------------------------------------------
    // Gene model loading (initialExtract3UTR)
    // --------------------------------------------------------------
    private static void loadGeneModel() throws IOException {
        Map<String, String> symbolMap = new HashMap<>();
        if (geneSymbolFile != null) {
            try (BufferedReader br = Files.newBufferedReader(Paths.get(geneSymbolFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    line = line.trim();
                    if (line.isEmpty() || line.startsWith("#")) continue;
                    String[] f = line.split("\t");
                    symbolMap.put(f[0], f[1]);
                }
            }
        }

        Map<String, List<UtrInfo>> chrMap = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(geneModelFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] f = line.split("\t");
                String chr = f[0];
                if (chr.contains("_")) continue;
                String refseq = f[3];
                String strand = f[5];
                String symbol = symbolMap.getOrDefault(refseq, "NA");
                String utrId = refseq + "|" + symbol + "|" + chr + "|" + strand;

                int txStart = Integer.parseInt(f[1]) + 1; // 1-based inclusive
                int txEnd = Integer.parseInt(f[2]);
                String[] relPosStr = f[f.length - 1].split(",");
                String[] lenStr = f[f.length - 2].split(",");
                int[] exonStarts = new int[relPosStr.length];
                int[] exonEnds = new int[relPosStr.length];
                for (int i = 0; i < relPosStr.length; i++) {
                    exonStarts[i] = Integer.parseInt(relPosStr[i]) + txStart;
                    exonEnds[i] = exonStarts[i] + Integer.parseInt(lenStr[i]) - 1;
                }

                int utrStart, utrEnd;
                if (strand.equals("+")) {
                    utrEnd = txEnd;
                    utrStart = txStart + Integer.parseInt(relPosStr[relPosStr.length - 1]) + 1;
                } else {
                    utrStart = Math.max(1, txStart);
                    // Perl: $fields[1] + ((split(/,/, $fields[10]))[0])   (likely blockSizes first element)
                    // We keep aligned with Perl's intended logic (which matches BED12 blockSizes[0])
                    utrEnd = txStart + Integer.parseInt(lenStr[0]);
                }

                UtrInfo utr = new UtrInfo(chr, utrStart, utrEnd, strand, exonStarts, exonEnds);
                utr.utrId = utrId;
                utrMap.put(utrId, utr);
                chrMap.computeIfAbsent(chr, k -> new ArrayList<>()).add(utr);
            }
        }

        // Pre-extend UTRs to avoid overlaps (Perl logic aligned)
        for (Map.Entry<String, List<UtrInfo>> e : chrMap.entrySet()) {
            List<UtrInfo> list = e.getValue();
            list.sort(Comparator.comparingInt(u -> u.start));
            int currLargest = 0;
            for (int i = 0; i < list.size(); i++) {
                UtrInfo u = list.get(i);
                if (u.strand.equals("+")) {
                    if (u.end < currLargest) continue;
                    currLargest = u.end;
                    if (i == list.size() - 1) {
                        u.end = u.annotatedEnd + extendSize;
                    } else {
                        // Find next UTR with end > current end
                        UtrInfo next = null;
                        for (int j = i + 1; j < list.size(); j++) {
                            UtrInfo candidate = list.get(j);
                            if (candidate.end > u.end) {
                                next = candidate;
                                break;
                            }
                        }
                        if (next == null) {
                            u.end = u.annotatedEnd + extendSize;
                        } else {
                            if (u.end < next.start) {
                                if (u.end + extendSize >= next.start) {
                                    u.end = next.start - 1;
                                } else {
                                    u.end = u.annotatedEnd + extendSize;
                                }
                            } else {
                                // u.end >= next.start : overlap already, keep current end (Perl: change_tick=1, last)
                                // Do nothing, u.end remains as is.
                                boolean ignored = false;
                            }
                        }
                    }
                } else {
                    if (i == 0) {
                        u.start = Math.max(1, u.annotatedStart - extendSize);
                        currLargest = u.end;
                    } else {
                        UtrInfo prev = list.get(i - 1);
                        int maxPrev = Math.max(prev.end, currLargest);
                        if (u.start > maxPrev) {
                            if (u.start - extendSize <= maxPrev) {
                                u.start = maxPrev + 1;
                            } else {
                                u.start = u.annotatedStart - extendSize;
                            }
                        }
                        if (u.end > currLargest) currLargest = u.end;
                    }
                }
            }
        }

        // Populate chrUtrMap for overlap removal later
        for (UtrInfo u : utrMap.values()) {
            chrUtrMap.computeIfAbsent(u.chr, k -> new LinkedHashMap<>())
                    .put(u.utrId, u);
        }
    }

    // --------------------------------------------------------------
    // Bedgraph loading
    // --------------------------------------------------------------
    private static void loadBedgraphFiles() throws Exception {
        ExecutorService exec = Executors.newFixedThreadPool(Math.min(inputFiles.size(), numThreads));
        List<Future<SampleData>> futures = new ArrayList<>();
        for (String f : inputFiles) {
            futures.add(exec.submit(() -> loadSingleBedgraph(f)));
        }
        for (Future<SampleData> f : futures) samples.add(f.get());
        exec.shutdown();

        System.err.println("Extracting UTR coverages...");
        for (SampleData sample : samples) {
            for (UtrInfo utr : utrMap.values()) {
                CoverageArray cov = sample.chrCoverage.get(utr.chr);
                double[] bpCov = (cov != null) ? extractUtrCoverage(utr, cov) : new double[utr.end - utr.start + 1];
                utr.sampleCoverages.add(bpCov);
            }
        }
    }

    private static SampleData loadSingleBedgraph(String path) throws IOException {
        SampleData data = new SampleData();
        double totalDepth = 0;
        Map<String, List<Integer>> startsMap = new HashMap<>();
        Map<String, List<Double>> covsMap = new HashMap<>();

        try (BufferedReader br = Files.newBufferedReader(Paths.get(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] f = line.split("\t");
                String chr = f[0];
                int start = Integer.parseInt(f[1]);
                int end = Integer.parseInt(f[2]);
                double cov = Double.parseDouble(f[f.length - 1]);
                totalDepth += cov * (end - start);

                List<Integer> st = startsMap.computeIfAbsent(chr, k -> new ArrayList<>());
                List<Double> cv = covsMap.computeIfAbsent(chr, k -> new ArrayList<>());
                if (!st.isEmpty() && start > st.get(st.size() - 1)) {
                    st.add(start);
                    cv.add(0.0);
                }
                st.add(end);
                cv.add(cov);
            }
        }

        for (String chr : startsMap.keySet()) {
            List<Integer> st = startsMap.get(chr);
            List<Double> cv = covsMap.get(chr);
            cv.add(0.0);
            int[] sArr = st.stream().mapToInt(Integer::intValue).toArray();
            double[] cArr = new double[cv.size()];
            for (int i = 0; i < cArr.length; i++) cArr[i] = cv.get(i);
            data.chrCoverage.put(chr, new CoverageArray(sArr, cArr));
        }
        data.totalDepth = totalDepth;
        return data;
    }

    private static double[] extractUtrCoverage(UtrInfo utr, CoverageArray cov) {
        int len = utr.end - utr.start + 1;
        double[] bp = new double[len];
        int[] starts = cov.starts;
        double[] covs = cov.covs;

        int leftIdx = Arrays.binarySearch(starts, utr.start);
        if (leftIdx < 0) leftIdx = -leftIdx - 2;
        leftIdx = Math.max(0, leftIdx);
        int rightIdx = Arrays.binarySearch(starts, utr.end);
        if (rightIdx < 0) rightIdx = -rightIdx - 2;
        rightIdx = Math.max(leftIdx, rightIdx);

        int rel = utr.start;
        for (int i = leftIdx; i <= rightIdx; i++) {
            int s = Math.max(utr.start, starts[i]) - rel;
            int e = Math.min(utr.end, (i < starts.length - 1 ? starts[i + 1] : utr.end)) - rel;
            double v = covs[i];
            for (int p = s; p < e; p++) bp[p] = v;
        }
        if (utr.strand.equals("-")) {
            for (int i = 0; i < len / 2; i++) {
                double t = bp[i];
                bp[i] = bp[len - 1 - i];
                bp[len - 1 - i] = t;
            }
        }
        return bp;
    }

    // --------------------------------------------------------------
    // Distal UTR identification (parallel)
    // --------------------------------------------------------------
    private static void identifyDistalUTRs() throws InterruptedException, ExecutionException {
        List<UtrInfo> utrList = new ArrayList<>(utrMap.values());
        int total = utrList.size();
        ExecutorService exec = Executors.newFixedThreadPool(numThreads);
        List<Future<?>> futures = new ArrayList<>();

        processedUtrs.set(0);
        done = false;
        startTime = System.currentTimeMillis();
        Thread progressThread = startProgressThread(total);

        for (UtrInfo utr : utrList) {
            futures.add(exec.submit(() -> {
                int distal = 0;
                for (double[] cov : utr.sampleCoverages) {
                    int d = findDistalEnd(utr, cov);
                    if (d > distal) distal = d;
                }
                if (distal > 0) {
                    if (utr.strand.equals("+")) {
                        utr.end = utr.start + distal - 1;
                    } else {
                        utr.start = utr.end - distal + 1;
                    }
                    utr.extended = true;
                } else {
                    utr.extended = false;
                }
                processedUtrs.incrementAndGet();
            }));
        }

        for (Future<?> f : futures) f.get();
        exec.shutdown();
        done = true;
        progressThread.join(1000);
    }

    private static Thread startProgressThread(int total) {
        Thread t = new Thread(() -> {
            while (!done) {
                try { Thread.sleep(1000); } catch (InterruptedException e) { break; }
                int p = processedUtrs.get();
                if (p == 0) continue;
                long elapsed = System.currentTimeMillis() - startTime;
                double rate = p * 1000.0 / elapsed;
                long remain = (long) ((total - p) / rate) * 1000;
                System.err.printf("\rProcessing UTRs: %d/%d (%.1f%%) ETA: %s",
                        p, total, p*100.0/total, formatDuration(remain));
            }
            System.err.printf("\rProcessing UTRs: %d/%d (100.0%%) Done.                    \n", total, total);
        });
        t.setDaemon(true);
        t.start();
        return t;
    }

    private static String formatDuration(long ms) {
        if (ms <= 0) return "0s";
        long s = ms / 1000;
        return String.format("%d:%02d:%02d", s/3600, (s%3600)/60, s%60);
    }

    // Perl-aligned distal end detection (using double precision)
    private static int findDistalEnd(UtrInfo utr, double[] cov) {
        int len = cov.length;
        if (len <= 1) return 0;

        int firstWin = Math.min(100, len);
        double firstMean = 0;
        for (int i = 0; i < firstWin; i++) firstMean += cov[i];
        firstMean /= firstWin;
        if (firstMean < 10) return 0;

        double[] prefix = new double[len + 1];
        for (int i = 0; i < len; i++) prefix[i + 1] = prefix[i] + cov[i];

        int annotatedLen = utr.annotatedEnd - utr.annotatedStart + 1;
        double largestMean = firstMean;
        int newEnd = 0;

        // First pass
        for (int i = 1; i < len; i++) {
            double frontMean;
            if (i >= windowSize) {
                int frontStart = i - windowSize;
                double frontSum = prefix[i] - prefix[frontStart];
                frontMean = frontSum / windowSize;
            } else {
                frontMean = firstMean;
            }
            if (frontMean > largestMean) largestMean = frontMean;

            int back1End = Math.min(i + windowSize, len);
            double back1Sum = prefix[back1End] - prefix[i];
            double back1Mean = back1Sum / (back1End - i);

            double back2Mean = 0;
            if (i + windowSize < len) {
                int back2End = Math.min(i + 2 * windowSize, len);
                double back2Sum = prefix[back2End] - prefix[i + windowSize];
                back2Mean = back2Sum / (back2End - i - windowSize);
            }

            double currCutoff = Math.max(largestMean * coverageCutoff, 1.0);
            // Align with Perl: threshold applies from the first base after annotated region
            double currPercentCutoff = (i >= annotatedLen) ? percentageCutoff : 0.0;

            int invalid1 = 0, total1 = back1End - i;
            for (int p = i; p < back1End; p++) if (cov[p] < currCutoff) invalid1++;
            int total2 = (i + windowSize < len) ? Math.min(i + 2 * windowSize, len) - (i + windowSize) : 0;
            int totalCombined = total1 + total2;

            boolean cond1 = (back2Mean < back1Mean || back2Mean < currCutoff) && back1Mean < currCutoff;
            boolean cond2 = invalid1 > (1 - currPercentCutoff) * totalCombined;

            if ((cond1 || cond2) && cov[i] < currCutoff) {
                return i;
            } else {
                if (back1Mean > frontMean * 4 && back1Mean > 10 && cov[i] > cov[i - 1] * 4 && i >= annotatedLen) {
                    newEnd = i - 1;
                    break;
                }
            }
        }

        if (newEnd == 0) return len;

        // Rescan up to newEnd
        largestMean = firstMean;
        for (int i = 1; i <= newEnd; i++) {
            double frontMean;
            if (i >= windowSize) {
                int frontStart = i - windowSize;
                double frontSum = prefix[i] - prefix[frontStart];
                frontMean = frontSum / windowSize;
            } else {
                frontMean = firstMean;
            }
            if (frontMean > largestMean) largestMean = frontMean;

            int back1End = Math.min(i + windowSize, newEnd + 1);
            double back1Sum = prefix[back1End] - prefix[i];
            double back1Mean = back1Sum / (back1End - i);

            double back2Mean = 0;
            if (i + windowSize <= newEnd) {
                int back2End = Math.min(i + 2 * windowSize, newEnd + 1);
                double back2Sum = prefix[back2End] - prefix[i + windowSize];
                back2Mean = back2Sum / (back2End - i - windowSize);
            }

            double currCutoff = Math.max(largestMean * coverageCutoff, 1.0);
            double currPercentCutoff = (i >= annotatedLen) ? percentageCutoff : 0.0;

            int invalid1 = 0, total1 = back1End - i;
            for (int p = i; p < back1End; p++) if (cov[p] < currCutoff) invalid1++;
            int total2 = (i + windowSize <= newEnd) ? Math.min(i + 2 * windowSize, newEnd + 1) - (i + windowSize) : 0;
            int totalCombined = total1 + total2;

            boolean cond1 = (back2Mean < back1Mean || back2Mean < currCutoff) && back1Mean < currCutoff;
            boolean cond2 = invalid1 > (1 - currPercentCutoff) * totalCombined;

            if ((cond1 || cond2) && cov[i] < currCutoff) {
                return i;
            }
        }
        return newEnd;
    }

    // --------------------------------------------------------------
    // Overlap removal (exactly as Perl)
    // --------------------------------------------------------------
    private static void removeOverlaps() {
        for (String chr : chrUtrMap.keySet()) {
            Map<String, UtrInfo> map = chrUtrMap.get(chr);
            List<String> ids = new ArrayList<>(map.keySet());
            ids.sort((a, b) -> {
                UtrInfo ua = map.get(a), ub = map.get(b);
                int cmp = Integer.compare(ua.start, ub.start);
                return cmp != 0 ? cmp : Integer.compare(ua.end, ub.end);
            });

            for (int i = 0; i < ids.size(); i++) {
                UtrInfo ui = map.get(ids.get(i));
                if (!ui.extended) continue;

                // Check previous UTRs
                for (int f = i - 1; f >= 0; f--) {
                    UtrInfo uf = map.get(ids.get(f));
                    if (ui.start > uf.end) {
                        if (f > i - 5) continue;
                        break;
                    }

                    // Perl: check if completely outside exons
                    boolean outside = true;
                    for (int exEnd : uf.exonEnds) if (exEnd >= ui.start) { outside = false; break; }
                    for (int exStart : uf.exonStarts) if (exStart <= ui.end) { outside = false; break; }
                    if (outside) continue;

                    // Overlap detected
                    if (ui.strand.equals(uf.strand)) {
                        if (ui.strand.equals("+") && uf.exonEnds.length > 1 && ui.start > uf.exonEnds[uf.exonEnds.length - 2])
                            continue;
                        if (ui.strand.equals("-") && uf.exonStarts.length > 1 && ui.end < uf.exonStarts[1])
                            continue;
                    } else {
                        if (ui.strand.equals("+") && uf.exonStarts[0] > ui.annotatedEnd &&
                                (ui.end - uf.exonStarts[0]) / (double) (ui.end - ui.start + 1) < 0.2) {
                            ui.end = uf.exonStarts[0] - 1;
                            continue;
                        }
                        if (ui.strand.equals("-") && uf.exonEnds[uf.exonEnds.length - 1] < ui.annotatedStart &&
                                (uf.exonEnds[uf.exonEnds.length - 1] - ui.start) / (double) (ui.end - ui.start + 1) < 0.2) {
                            ui.start = uf.exonEnds[uf.exonEnds.length - 1] + 1;
                            continue;
                        }
                    }
                    ui.extended = false;
                    break;
                }
                if (!ui.extended) continue;

                // Check next UTRs
                for (int f = i + 1; f < ids.size(); f++) {
                    UtrInfo uf = map.get(ids.get(f));
                    if (ui.end < uf.start) {
                        if (f < i + 5) continue;
                        break;
                    }

                    // Perl: check if completely outside exons
                    boolean outside = true;
                    for (int exEnd : uf.exonEnds) if (exEnd >= ui.start) { outside = false; break; }
                    for (int exStart : uf.exonStarts) if (exStart <= ui.end) { outside = false; break; }
                    if (outside) continue;

                    if (ui.strand.equals(uf.strand)) {
                        if (ui.strand.equals("+") && uf.exonEnds.length > 1 && ui.start > uf.exonEnds[uf.exonEnds.length - 2])
                            continue;
                        if (ui.strand.equals("-") && uf.exonStarts.length > 1 && ui.end < uf.exonStarts[1])
                            continue;
                    } else {
                        if (ui.strand.equals("+") && uf.exonStarts[0] > ui.annotatedEnd &&
                                (ui.end - uf.exonStarts[0]) / (double) (ui.end - ui.start + 1) < 0.2) {
                            ui.end = uf.exonStarts[0] - 1;
                            continue;
                        }
                        if (ui.strand.equals("-") && uf.exonEnds[uf.exonEnds.length - 1] < ui.annotatedStart &&
                                (uf.exonEnds[uf.exonEnds.length - 1] - ui.start) / (double) (ui.end - ui.start + 1) < 0.2) {
                            ui.start = uf.exonEnds[uf.exonEnds.length - 1] + 1;
                            continue;
                        }
                    }
                    ui.extended = false;
                    break;
                }
            }
        }
    }

    private static void writeOutput() throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
            for (UtrInfo u : utrMap.values()) {
                if (!u.extended) continue;
                pw.printf("%s\t%d\t%d\t%s\t0\t%s\n",
                        u.chr, u.start, u.end, u.utrId, u.strand);
            }
        }
    }

    // --------------------------------------------------------------
    // Data classes (double precision)
    // --------------------------------------------------------------
    static class UtrInfo {
        String chr, strand, utrId;
        int start, end;
        int annotatedStart, annotatedEnd;
        int[] exonStarts, exonEnds;
        boolean extended = true;
        List<double[]> sampleCoverages = new ArrayList<>();

        UtrInfo(String chr, int s, int e, String strand, int[] exonStarts, int[] exonEnds) {
            this.chr = chr;
            this.annotatedStart = s;
            this.annotatedEnd = e;
            this.start = s;
            this.end = e;
            this.strand = strand;
            this.exonStarts = exonStarts;
            this.exonEnds = exonEnds;
        }
    }

    static class SampleData {
        double totalDepth;
        Map<String, CoverageArray> chrCoverage = new HashMap<>();
    }

    static class CoverageArray {
        int[] starts;
        double[] covs;
        CoverageArray(int[] s, double[] c) { starts = s; covs = c; }
    }
}