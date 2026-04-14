package com.vbuser.apa;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * High-performance Java 8 reimplementation of predictAPA.pl.
 * Strictly algorithm-equivalent to Perl version with all bugs preserved.
 * Optimized for CPU efficiency: prefix sums, array views, minimal allocation.
 * Bounds-checked for robustness.
 * Fixed mathematical discrepancies compared to Perl:
 * 1. sample_marks skipping logic implemented (invalid samples excluded from MSE/abundance)
 * 2. Removed Math.max(1, ...) for window/minDist calculations to match Perl int division
 * 3. Negative strand breakpoints no longer sorted (preserve discovery order)
 * 4. Mean denominator uses actual slice length (already correct)
 * MODIFICATIONS:
 * - Commented out coverageThreshold sample invalidation to match Perl behavior (reduces NAs)
 * - Upgraded all coverage values from float to double for higher precision
 */
@SuppressWarnings({"CallToPrintStackTrace", "UseOfSystemOutOrSystemErr"})
public class PredictAPA {

    // Command line parameters
    private static List<String> inputFiles;
    private static int numGroups;
    private static int[] groupNums;
    private static String utrFile;
    private static String outputFile;
    private static double dropCutoff = 0.2;
    private static double coverageThreshold = 20.0;   // NOT USED due to Perl bug
    private static int apaMinDist = 100;
    private static int scanWindowSize = 50;

    // Data structures
    private static Map<String, UtrInfo> utrMap = new LinkedHashMap<>();
    private static final List<SampleData> samples = new ArrayList<>();

    // Statistics
    private static final AtomicInteger writtenCount = new AtomicInteger(0);
    private static final AtomicInteger noApaSiteCount = new AtomicInteger(0);
    private static final AtomicInteger zeroTotalExpCount = new AtomicInteger(0);
    private static final AtomicInteger shortUtrCount = new AtomicInteger(0);

    // Progress
    private static final AtomicInteger processedCount = new AtomicInteger(0);
    private static volatile boolean done = false;
    private static long startTime;

    public static void main(String[] args) throws Exception {
        parseArguments(args);
        validateArguments();

        System.err.println("Loading UTR annotation...");
        loadUtrAnnotation();

        System.err.println("Loading bedgraph files in parallel...");
        loadBedgraphFiles();

        double meanDepth = samples.stream().mapToDouble(s -> s.totalDepth).average().orElse(1.0);
        double[] sampleWeights = samples.stream().mapToDouble(s -> s.totalDepth / meanDepth).toArray();

        System.err.println("Starting de novo APA estimation...");

        List<Map.Entry<String, UtrInfo>> validUtrs = utrMap.entrySet().stream()
                .filter(e -> e.getValue().sampleIntervals.size() == samples.size())
                .collect(Collectors.toList());

        int totalUtrs = validUtrs.size();

        // Write header
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            List<String> header = new ArrayList<>();
            header.add("Gene");
            header.add("Mean_Squared_Error");
            header.add("Predicted_APA");
            header.add("Loci");
            for (int g = 0; g < numGroups; g++) {
                for (int rep = 0; rep < groupNums[g]; rep++) {
                    header.add("Group_" + (g + 1) + "_" + (rep + 1) + "_Separate_Exp");
                    header.add("Group_" + (g + 1) + "_" + (rep + 1) + "_Total_Exp");
                }
            }
            writer.println(String.join("\t", header));
        }

        int threads = Math.min(Runtime.getRuntime().availableProcessors(), 40);
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        List<Future<String>> futures = new ArrayList<>();

        BlockingQueue<String> resultQueue = new LinkedBlockingQueue<>();
        Thread writerThread = new Thread(() -> {
            try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile, true))) {
                while (true) {
                    String line = resultQueue.poll(100, TimeUnit.MILLISECONDS);
                    if (line == null) {
                        if (futures.stream().allMatch(Future::isDone) && resultQueue.isEmpty()) break;
                        continue;
                    }
                    writer.println(line);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        writerThread.start();

        processedCount.set(0);
        done = false;
        startTime = System.currentTimeMillis();
        Thread progressThread = getThread(totalUtrs);

        for (Map.Entry<String, UtrInfo> entry : validUtrs) {
            futures.add(executor.submit(() -> processUtr(entry, sampleWeights, resultQueue)));
        }

        for (Future<String> future : futures) {
            try { future.get(); } catch (Exception e) { e.printStackTrace(); }
        }

        done = true;
        progressThread.interrupt();
        progressThread.join(1000);
        executor.shutdown();
        writerThread.join();

        System.err.println("\n--- Skip Statistics ---");
        System.err.println("  written: " + writtenCount.get());
        System.err.println("  no_apa_site: " + noApaSiteCount.get());
        System.err.println("  zero_total_exp: " + zeroTotalExpCount.get());
        System.err.println("  too_short: " + shortUtrCount.get());
        System.err.println("Done!");
    }

    private static Thread getThread(int totalUtrs) {
        Thread progressThread = new Thread(() -> {
            while (!done) {
                try { Thread.sleep(1000); } catch (InterruptedException e) { break; }
                int processed = processedCount.get();
                if (processed == 0) continue;
                long elapsed = System.currentTimeMillis() - startTime;
                double rate = processed * 1000.0 / elapsed;
                long remaining = (long) ((totalUtrs - processed) / rate) * 1000;
                String eta = formatDuration(remaining);
                double percent = processed * 100.0 / totalUtrs;
                System.err.printf("\rProcessing UTRs: %d/%d (%.1f%%) | ETA: %s",
                        processed, totalUtrs, percent, eta);
            }
            System.err.printf("\rProcessing UTRs: %d/%d (100.0%%) | Done.                    \n",
                    totalUtrs, totalUtrs);
        });
        progressThread.setDaemon(true);
        progressThread.start();
        return progressThread;
    }

    private static String formatDuration(long millis) {
        if (millis <= 0) return "0s";
        long seconds = millis / 1000;
        long minutes = seconds / 60;
        seconds %= 60;
        long hours = minutes / 60;
        minutes %= 60;
        if (hours > 0) return String.format("%dh%02dm", hours, minutes);
        else if (minutes > 0) return String.format("%dm%02ds", minutes, seconds);
        else return String.format("%ds", seconds);
    }

    private static void parseArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-i":
                    inputFiles = new ArrayList<>();
                    while (i + 1 < args.length && !args[i + 1].startsWith("-")) inputFiles.add(args[++i]);
                    break;
                case "-g": numGroups = Integer.parseInt(args[++i]); break;
                case "-n":
                    List<Integer> nums = new ArrayList<>();
                    while (i + 1 < args.length && !args[i + 1].startsWith("-")) nums.add(Integer.parseInt(args[++i]));
                    groupNums = nums.stream().mapToInt(Integer::intValue).toArray();
                    break;
                case "-u": utrFile = args[++i]; break;
                case "-o": outputFile = args[++i]; break;
                case "-d": dropCutoff = Double.parseDouble(args[++i]); break;
                case "-c": coverageThreshold = Double.parseDouble(args[++i]); break;
                case "-a": apaMinDist = Integer.parseInt(args[++i]); break;
                case "-w": scanWindowSize = Integer.parseInt(args[++i]); break;
            }
        }
    }

    private static void validateArguments() {
        if (inputFiles == null || inputFiles.isEmpty()) error("-i required");
        if (numGroups <= 0) error("-g must be >= 1");
        if (groupNums == null || groupNums.length != numGroups) error("-n must have " + numGroups + " values");
        if (inputFiles.size() != Arrays.stream(groupNums).sum()) error("number of input files != sum of -n");
        if (utrFile == null) error("-u required");
        if (outputFile == null) error("-o required");
        if (dropCutoff <= 0 || dropCutoff >= 1) error("-d in (0,1)");
        if (coverageThreshold < 10) error("-c >=10");
        if (apaMinDist < 20) error("-a >=20");
        if (scanWindowSize < 20) error("-w >=20");
    }

    private static void error(String msg) {
        throw new IllegalArgumentException("PredictAPA error: " + msg);
    }

    private static void loadUtrAnnotation() throws IOException {
        try (BufferedReader br = Files.newBufferedReader(Paths.get(utrFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] parts = line.split("\t");
                String chr = parts[0];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                String strand = parts[parts.length - 1];
                if (start + 50 < end) {
                    String gene = parts[3];
                    utrMap.put(gene, new UtrInfo(chr, start, end, strand, chr + ":" + start + "-" + end));
                }
            }
        }
        utrMap = utrMap.entrySet().stream()
                .sorted((a, b) -> {
                    UtrInfo u1 = a.getValue(), u2 = b.getValue();
                    int cmp = u1.chr.compareTo(u2.chr);
                    return cmp != 0 ? cmp : Integer.compare(u1.start, u2.start);
                })
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
    }

    private static void loadBedgraphFiles() throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(Math.min(inputFiles.size(),
                Math.min(Runtime.getRuntime().availableProcessors(), 40)));
        List<Future<SampleData>> futures = new ArrayList<>();
        for (String f : inputFiles) futures.add(executor.submit(() -> loadSingleBedgraph(f)));
        for (Future<SampleData> f : futures) samples.add(f.get());
        executor.shutdown();

        System.err.println("Extracting UTR coverages...");
        for (SampleData sample : samples) {
            Map<String, CoverageArray> covMap = sample.chrCoverage;
            String preChr = "";
            CoverageArray currCov = null;
            int leftOffset = 0;
            for (UtrInfo utr : utrMap.values()) {
                String chr = utr.chr;
                if (!covMap.containsKey(chr)) {
                    utr.sampleIntervals.add(new UtrInterval(utr.start, utr.end, utr.strand, new int[]{utr.start}, new double[]{0.0}));
                    continue;
                }
                if (!chr.equals(preChr)) {
                    currCov = covMap.get(chr);
                    preChr = chr;
                } else if (leftOffset > 0) {
                    currCov = currCov.slice(leftOffset);
                }
                assert currCov != null;
                int[] starts = currCov.starts;
                double[] covs = currCov.covs;

                int leftIdx = Arrays.binarySearch(starts, utr.start);
                if (leftIdx < 0) leftIdx = -leftIdx - 2;
                leftIdx = Math.max(0, leftIdx);
                int rightIdx = Arrays.binarySearch(starts, utr.end);
                if (rightIdx < 0) rightIdx = -rightIdx - 2;
                rightIdx = Math.max(leftIdx, rightIdx);

                int len = rightIdx - leftIdx + 1;
                int[] extStarts = new int[len];
                double[] extCovs = new double[len];
                System.arraycopy(starts, leftIdx, extStarts, 0, len);
                System.arraycopy(covs, leftIdx, extCovs, 0, len);

                if (extStarts[0] > utr.start) {
                    extStarts = prependInt(extStarts, utr.start);
                    extCovs = prependDouble(extCovs);
                }
                if (extStarts[extStarts.length - 1] < utr.end) {
                    extStarts = appendInt(extStarts, utr.end);
                    extCovs = appendDouble(extCovs, extCovs[extCovs.length - 1]);
                } else {
                    extStarts[extStarts.length - 1] = utr.end;
                }
                utr.sampleIntervals.add(new UtrInterval(utr.start, utr.end, utr.strand, extStarts, extCovs));
                leftOffset = leftIdx;
            }
        }
    }

    private static SampleData loadSingleBedgraph(String path) throws IOException {
        SampleData data = new SampleData();
        double totalDepth = 0.0;
        Map<String, List<Integer>> startsMap = new HashMap<>();
        Map<String, List<Double>> covsMap = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] parts = line.split("\t");
                String chr = parts[0];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                double cov = Double.parseDouble(parts[parts.length - 1]);
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

    private static int[] prependInt(int[] a, int v) { int[] r = new int[a.length+1]; r[0]=v; System.arraycopy(a,0,r,1,a.length); return r; }
    private static double[] prependDouble(double[] a) { double[] r = new double[a.length+1]; r[0]= 0.0; System.arraycopy(a,0,r,1,a.length); return r; }
    private static int[] appendInt(int[] a, int v) { int[] r = Arrays.copyOf(a, a.length+1); r[a.length]=v; return r; }
    private static double[] appendDouble(double[] a, double v) { double[] r = Arrays.copyOf(a, a.length+1); r[a.length]=v; return r; }

    // --------------------------------------------------------------
    // Core optimized processing
    // --------------------------------------------------------------
    private static String processUtr(Map.Entry<String, UtrInfo> entry, double[] weights, BlockingQueue<String> queue) {
        UtrInfo utr = entry.getValue();
        int len = utr.end - utr.start + 1;
        if (len < 100) {
            shortUtrCount.incrementAndGet();
            processedCount.incrementAndGet();
            return null;
        }

        int nSamples = samples.size();
        double[][] rawCovs = new double[nSamples][];
        for (int s = 0; s < nSamples; s++) {
            rawCovs[s] = utr.sampleIntervals.get(s).toBpCoverage();
            if (rawCovs[s].length != len) {
                // In case of unexpected length mismatch, pad/truncate (shouldn't happen)
                if (rawCovs[s].length < len) {
                    double[] tmp = new double[len];
                    System.arraycopy(rawCovs[s], 0, tmp, 0, rawCovs[s].length);
                    rawCovs[s] = tmp;
                } else if (rawCovs[s].length > len) {
                    rawCovs[s] = Arrays.copyOf(rawCovs[s], len);
                }
            }
        }

        // Identify internal APA regions and sample validity marks
        RegionsResult regResult = identifyRegions(rawCovs, len);
        if (regResult.regions.isEmpty()) {
            noApaSiteCount.incrementAndGet();
            processedCount.incrementAndGet();
            return null;
        }

        boolean[] validSample = regResult.sampleMarks;

        // Normalize coverage by sample weights
        double[][] normCovs = new double[nSamples][];
        for (int s = 0; s < nSamples; s++) {
            double[] raw = rawCovs[s];
            double[] norm = new double[len];
            double w = weights[s];
            for (int i = 0; i < len; i++) norm[i] = raw[i] / w;
            normCovs[s] = norm;
        }

        Candidate best = findBestBreakpoints(regResult.regions, normCovs, len, nSamples, validSample);
        if (best == null) {
            noApaSiteCount.incrementAndGet();
            processedCount.incrementAndGet();
            return null;
        }

        // Check that total expression > 0 for all valid samples
        for (int s = 0; s < nSamples; s++) {
            if (!validSample[s]) continue;
            double total = 0.0;
            for (double v : best.abundances[s]) total += v;
            if (total <= 0) {
                zeroTotalExpCount.incrementAndGet();
                processedCount.incrementAndGet();
                return null;
            }
        }

        // Convert relative breakpoints to genomic coordinates (preserve discovery order)
        int[] genomicBreaks = new int[best.breakPoints.length];
        if (utr.strand.equals("+")) {
            for (int i = 0; i < genomicBreaks.length; i++)
                genomicBreaks[i] = best.breakPoints[i] + utr.start - 1;
        } else {
            for (int i = 0; i < genomicBreaks.length; i++)
                genomicBreaks[i] = utr.end - best.breakPoints[i] + 1;
        }

        StringBuilder sb = new StringBuilder();
        sb.append(entry.getKey()).append('\t');
        sb.append(String.format("%.1f", best.meanSquaredError)).append('\t');
        sb.append(Arrays.stream(genomicBreaks).mapToObj(Integer::toString).collect(Collectors.joining(","))).append('\t');
        sb.append(utr.loci);
        for (int s = 0; s < nSamples; s++) {
            if (!validSample[s]) {
                sb.append("\tNA\tNA");
            } else {
                double[] abund = best.abundances[s];
                String sep = Arrays.stream(abund).mapToObj(v -> String.format("%.2f", v)).collect(Collectors.joining(","));
                double total = Arrays.stream(abund).sum();
                sb.append('\t').append(sep).append('\t').append(String.format("%.2f", total));
            }
        }

        writtenCount.incrementAndGet();
        processedCount.incrementAndGet();
        queue.add(sb.toString());
        return sb.toString();
    }

    /**
     * Identify candidate internal APA regions, also marks samples as invalid
     * if they meet Perl's sample_marks clearing conditions.
     * Fix: removed coverageThreshold invalidation to match Perl behavior.
     */
    private static RegionsResult identifyRegions(double[][] rawCovs, int len) {
        double cutoff = 1.0 - dropCutoff;
        int window = (len < 4 * scanWindowSize) ? scanWindowSize / 2 : scanWindowSize;
        int minDist = (len < 3 * apaMinDist) ? apaMinDist / 2 : apaMinDist;
        int searchEnd = len - minDist;
        if (window >= searchEnd) {
            boolean[] marks = new boolean[rawCovs.length];
            Arrays.fill(marks, true);
            return new RegionsResult(Collections.emptyList(), marks);
        }

        int nSamples = rawCovs.length;
        boolean[] sampleMarks = new boolean[nSamples];
        Arrays.fill(sampleMarks, true);
        List<List<int[]>> allRegions = new ArrayList<>(nSamples);

        for (int s = 0; s < nSamples; s++) {
            double[] cov = rawCovs[s];
            double[] prefix = new double[len + 1];
            for (int i = 0; i < len; i++) prefix[i + 1] = prefix[i] + cov[i];

            int topLen = Math.min(100, len);
            double topMean = prefix[topLen] / topLen;

            List<int[]> regions = new ArrayList<>();
            int curr = window;
            boolean sampleValid = true;

            while (curr <= searchEnd) {
                // Front window mean
                int fwStart = Math.max(0, curr - window);
                int fwEnd = curr;
                int fwSize = fwEnd - fwStart;
                double frontWinMean = (fwSize > 0) ? ((prefix[fwEnd] - prefix[fwStart]) / fwSize) : 0.0;

                // Back window mean
                int bwStart = curr;
                int bwEnd = Math.min(curr + window, len);
                int bwSize = bwEnd - bwStart;
                double backWinMean = (bwSize > 0) ? ((prefix[bwEnd] - prefix[bwStart]) / bwSize) : 0.0;

                // Front mean up to curr
                double frontMean;
                if (regions.isEmpty()) {
                    frontMean = (curr > 0) ? (prefix[curr] / curr) : 0.0;
                } else {
                    int[] last = regions.get(regions.size() - 1);
                    int mid = (last[0] + last[1]) / 2;
                    mid = Math.max(0, Math.min(mid, curr));
                    int span = curr - mid;
                    frontMean = (span > 0) ? ((prefix[curr] - prefix[mid]) / span) : 0.0;
                }

                if (curr < Math.max(window, 0.2 * len) && backWinMean > topMean) {
                    topMean = backWinMean;
                }

                // Perl sample_marks invalidation condition
                if (backWinMean > topMean && frontMean * 4 < backWinMean) {
                    regions.clear();
                    sampleValid = false;
                    break;
                }

                double currCutoff = Math.max(frontMean, frontWinMean) * cutoff;
                double currCovValue = (curr < len) ? cov[curr] : 0.0;
                if (backWinMean < currCutoff && currCovValue < currCutoff &&
                        Math.max(frontMean, frontWinMean) - backWinMean >= 10) {

                    int remainingLen = len - curr;
                    int validBack = 0;
                    for (int i = curr; i < len; i++) if (cov[i] < currCutoff) validBack++;
                    boolean backOk = (remainingLen > 0) && (validBack / (double) remainingLen >= 0.75);
                    boolean allLess = true;
                    double maxFront = Math.max(frontMean, frontWinMean);
                    for (int i = curr; i < len; i++) {
                        if (cov[i] >= maxFront) {
                            allLess = false;
                            break;
                        }
                    }
                    if (backOk && allLess && remainingLen > 0.05 * len) {
                        int regStart = curr - minDist / 2;
                        int regEnd = curr + minDist / 2;
                        regions.add(new int[]{regStart, regEnd});
                        curr += window;
                        continue;
                    }
                }
                curr++;
            }

            sampleMarks[s] = sampleValid;
            // ** FIX: Commented out to match Perl behavior (prevents unnecessary NA) **
            // if (sampleValid && topMean < coverageThreshold) {
            //     sampleMarks[s] = false;
            // }
            allRegions.add(sampleValid ? regions : new ArrayList<>());
        }

        // Merge overlapping regions across samples
        List<int[]> merged = new ArrayList<>();
        for (List<int[]> regs : allRegions) merged.addAll(regs);
        if (merged.isEmpty()) return new RegionsResult(Collections.emptyList(), sampleMarks);
        merged.sort(Comparator.comparingInt(a -> a[0]));

        List<int[]> finalRegions = new ArrayList<>();
        int[] prev = null;
        int mergeDist = (len < 500) ? 0 : window;
        for (int[] r : merged) {
            if (prev == null) prev = r.clone();
            else if (prev[1] + mergeDist < r[0]) {
                finalRegions.add(prev);
                prev = r.clone();
            } else {
                prev[1] = Math.max(prev[1], r[1]);
            }
        }
        if (prev != null) finalRegions.add(prev);
        return new RegionsResult(finalRegions, sampleMarks);
    }

    private static Candidate findBestBreakpoints(List<int[]> regions, double[][] normCovs, int len, int nSamples, boolean[] validSample) {
        if (regions.size() > 1) return efficientSearch(regions, normCovs, len, nSamples, validSample);
        int[] region = regions.get(0);
        Candidate best = null;
        double bestMSE = Double.POSITIVE_INFINITY;
        for (int p = region[0]; p <= region[1]; p++) {
            int[] breaks = {p};
            double sumMSE = 0.0;
            List<double[]> abunds = new ArrayList<>();
            int validCount = 0;
            boolean ok = true;
            for (int s = 0; s < nSamples; s++) {
                if (!validSample[s]) continue;
                EstimateResult er = estimateAbundance(breaks, normCovs[s]);
                if (er == null) { ok = false; break; }
                sumMSE += er.mse;
                abunds.add(er.abund);
                validCount++;
            }
            if (!ok || validCount == 0) continue;
            double avg = sumMSE / validCount;
            if (avg < bestMSE) {
                bestMSE = avg;
                double[][] abundArray = new double[nSamples][];
                int idx = 0;
                for (int s = 0; s < nSamples; s++) {
                    if (validSample[s]) abundArray[s] = abunds.get(idx++);
                    else abundArray[s] = new double[0];
                }
                best = new Candidate(breaks, avg, abundArray);
            }
        }
        return best;
    }

    private static Candidate efficientSearch(List<int[]> regions, double[][] normCovs, int len, int nSamples, boolean[] validSample) {
        int nRegions = regions.size();
        int[] selected = new int[nRegions];
        for (int i = 0; i < nRegions; i++) {
            int[] reg = regions.get(i);
            int partStart, partEnd;
            if (i == 0) {
                partStart = 0;
                partEnd = (regions.get(i+1)[0] + regions.get(i+1)[1]) / 2;
            } else if (i == nRegions - 1) {
                partStart = selected[i-1];
                partEnd = len;
            } else {
                partStart = selected[i-1];
                partEnd = (regions.get(i+1)[0] + regions.get(i+1)[1]) / 2;
            }
            partStart = Math.max(0, Math.min(partStart, len));
            partEnd = Math.max(partStart, Math.min(partEnd, len));
            int relStart = reg[0] - partStart;
            int relEnd = reg[1] - partStart;
            relStart = Math.max(0, relStart);
            relEnd = Math.min(relEnd, partEnd - partStart);
            if (relStart > relEnd) continue;
            double bestMSE = Double.POSITIVE_INFINITY;
            int bestPoint = -1;
            for (int p = relStart; p <= relEnd; p++) {
                double sumMSE = 0.0;
                int validCount = 0;
                for (int s = 0; s < nSamples; s++) {
                    if (!validSample[s]) continue;
                    double[] cov = normCovs[s];
                    EstimateResult er = estimateAbundance(new int[]{p}, cov, partStart, partEnd);
                    if (er == null) { sumMSE = Double.POSITIVE_INFINITY; break; }
                    sumMSE += er.mse;
                    validCount++;
                }
                if (validCount == 0) continue;
                double avg = sumMSE / validCount;
                if (avg < bestMSE) {
                    bestMSE = avg;
                    bestPoint = p;
                }
            }
            if (bestPoint == -1) return null;
            selected[i] = bestPoint + partStart;
        }

        double sumMSE = 0.0;
        int validCount = 0;
        double[][] abundArray = new double[nSamples][];
        for (int s = 0; s < nSamples; s++) {
            if (!validSample[s]) continue;
            EstimateResult er = estimateAbundance(selected, normCovs[s]);
            if (er == null) return null;
            sumMSE += er.mse;
            abundArray[s] = er.abund;
            validCount++;
        }
        for (int s = 0; s < nSamples; s++) {
            if (!validSample[s]) abundArray[s] = new double[0];
        }
        if (validCount == 0) return null;
        return new Candidate(selected, sumMSE / validCount, abundArray);
    }

    private static EstimateResult estimateAbundance(int[] breaks, double[] cov, int start, int end) {
        int n = breaks.length;
        int[] segStarts = new int[n+1];
        int[] segEnds = new int[n+1];
        segStarts[0] = 0;
        segEnds[0] = breaks[0];
        for (int i = 1; i < n; i++) {
            segStarts[i] = breaks[i-1];
            segEnds[i] = breaks[i];
        }
        segStarts[n] = breaks[n-1];
        segEnds[n] = end - start;

        double[] means = new double[n+1];
        for (int i = 0; i <= n; i++) {
            int s = segStarts[i];
            int e = segEnds[i];
            if (e > s) {
                double sum = 0.0;
                int maxIdx = Math.min(start + e, cov.length);
                for (int j = start + s; j < maxIdx; j++) sum += cov[j];
                means[i] = sum / (e - s);
            }
        }
        double[] abund = means.clone();
        for (int i = n-1; i >= 0; i--) {
            double down = 0.0;
            for (int j = i+1; j <= n; j++) down += abund[j];
            abund[i] = Math.max(0.0, abund[i] - down);
        }
        double mse = 0.0;
        int cnt = 0;
        for (int i = 0; i <= n; i++) {
            double cum = 0.0;
            for (int j = i; j <= n; j++) cum += abund[j];
            int maxP = Math.min(start + segEnds[i], cov.length);
            for (int p = start + segStarts[i]; p < maxP; p++) {
                double diff = cov[p] - cum;
                mse += diff * diff;
                cnt++;
            }
        }
        if (cnt == 0) return null;
        return new EstimateResult(mse / cnt, abund);
    }

    private static EstimateResult estimateAbundance(int[] breaks, double[] cov) {
        int n = breaks.length;
        int len = cov.length;
        int[] segStarts = new int[n+1];
        int[] segEnds = new int[n+1];
        segStarts[0] = 0;
        segEnds[0] = Math.min(breaks[0], len);
        for (int i = 1; i < n; i++) {
            segStarts[i] = Math.min(breaks[i-1], len);
            segEnds[i] = Math.min(breaks[i], len);
        }
        segStarts[n] = Math.min(breaks[n-1], len);
        segEnds[n] = len;

        double[] means = new double[n+1];
        for (int i = 0; i <= n; i++) {
            int s = segStarts[i];
            int e = segEnds[i];
            if (e > s) {
                double sum = 0.0;
                for (int j = s; j < e; j++) sum += cov[j];
                means[i] = sum / (e - s);
            }
        }
        double[] abund = means.clone();
        for (int i = n-1; i >= 0; i--) {
            double down = 0.0;
            for (int j = i+1; j <= n; j++) down += abund[j];
            abund[i] = Math.max(0.0, abund[i] - down);
        }
        double mse = 0.0;
        int cnt = 0;
        for (int i = 0; i <= n; i++) {
            double cum = 0.0;
            for (int j = i; j <= n; j++) cum += abund[j];
            for (int p = segStarts[i]; p < segEnds[i]; p++) {
                double diff = cov[p] - cum;
                mse += diff * diff;
                cnt++;
            }
        }
        if (cnt == 0) return null;
        return new EstimateResult(mse / cnt, abund);
    }

    // --------------------------------------------------------------
    // Data classes (updated to use double for coverage)
    // --------------------------------------------------------------
    static class UtrInfo {
        String chr; int start, end; String strand; String loci;
        List<UtrInterval> sampleIntervals = new ArrayList<>();
        UtrInfo(String c, int s, int e, String str, String loc) { chr=c; start=s; end=e; strand=str; loci=loc; }
    }

    static class UtrInterval {
        int utrStart, utrEnd; String strand; int[] starts; double[] covs;
        UtrInterval(int s, int e, String str, int[] st, double[] cv) { utrStart=s; utrEnd=e; strand=str; starts=st; covs=cv; }

        double[] toBpCoverage() {
            int len = utrEnd - utrStart + 1;
            double[] bp = new double[len];
            int rel = utrStart;
            for (int i = 0; i < starts.length - 1; i++) {
                int s = starts[i] - rel;
                int e = starts[i+1] - rel - 1;
                double v = covs[i];
                for (int p = s; p <= e; p++) if (p >= 0 && p < len) bp[p] = v;
            }
            int last = starts.length - 1;
            int s = starts[last] - rel;
            int e = len - 1;
            double v = covs[last];
            for (int p = s; p <= e; p++) if (p >= 0) bp[p] = v;
            if (strand.equals("-")) {
                for (int i = 0; i < len/2; i++) {
                    double t = bp[i]; bp[i] = bp[len-1-i]; bp[len-1-i] = t;
                }
            }
            return bp;
        }
    }

    static class SampleData { double totalDepth; Map<String, CoverageArray> chrCoverage = new HashMap<>(); }

    static class CoverageArray {
        int[] starts; double[] covs;
        CoverageArray(int[] s, double[] c) { starts = s; covs = c; }
        CoverageArray slice(int off) {
            int nl = starts.length - off;
            int[] ns = new int[nl];
            double[] nc = new double[nl];
            System.arraycopy(starts, off, ns, 0, nl);
            System.arraycopy(covs, off, nc, 0, nl);
            return new CoverageArray(ns, nc);
        }
    }

    static class EstimateResult { double mse; double[] abund; EstimateResult(double m, double[] a) { mse=m; abund=a; } }

    static class Candidate {
        int[] breakPoints; double meanSquaredError; double[][] abundances;
        Candidate(int[] b, double m, double[][] a) { breakPoints=b; meanSquaredError=m; abundances=a; }
    }

    static class RegionsResult {
        final List<int[]> regions;
        final boolean[] sampleMarks;
        RegionsResult(List<int[]> r, boolean[] m) { regions = r; sampleMarks = m; }
    }
}