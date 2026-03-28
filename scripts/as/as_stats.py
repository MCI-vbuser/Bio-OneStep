#!/usr/bin/env python3
"""
AS data statistics script
Usage: python as_stats.py <expressed_dir> <enhanced_dir>
Generates stats.db (SQLite) and index.html (dark theme English UI).
"""

import os
import sys
import argparse
import glob
import csv
import json
import sqlite3
import re
import numpy as np
from collections import defaultdict
import pandas as pd

# ------------------------------------------------------------
# Parse .ioe files: count events per type from file names
def parse_ioe_files(ioe_dir):
    event_type_count = defaultdict(int)
    total_events = 0
    ioe_files = glob.glob(os.path.join(ioe_dir, "*.ioe"))
    for ioe_file in ioe_files:
        basename = os.path.basename(ioe_file)
        match = re.search(r'[._]?events?_([A-Z0-9_]+)\.ioe', basename)
        if match:
            etype = match.group(1)
        else:
            etype = os.path.splitext(basename)[0]
            if etype.startswith('events_'):
                etype = etype[7:]
            if not etype:
                etype = "unknown"
        with open(ioe_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            try:
                header = next(reader)
            except StopIteration:
                continue
            row_count = 0
            for row in reader:
                if row and not row[0].startswith('#'):
                    row_count += 1
        event_type_count[etype] += row_count
        total_events += row_count
    return total_events, dict(event_type_count)

# ------------------------------------------------------------
# Parse .psi files with detailed quality metrics
def parse_psi_files(psi_dir):
    psi_files = glob.glob(os.path.join(psi_dir, "*.psi"))
    if not psi_files:
        return None

    all_psi_values = []
    event_mean_psis = []
    event_valid_counts = []
    type_stats_list = []
    total_events_all = 0
    total_cells_all = 0
    valid_cells_all = 0
    zero_valid_events = 0

    for psi_file in psi_files:
        basename = os.path.basename(psi_file)
        match = re.match(r'([A-Z0-9]+)\.psi', basename)
        etype = match.group(1) if match else "unknown"
        with open(psi_file, 'r') as f:
            lines = f.readlines()
        if not lines:
            continue
        events_in_file = 0
        valid_events_in_file = 0
        valid_cells_in_file = 0
        file_psi_values = []
        file_event_means = []
        file_event_valid_counts = []
        file_total_cells = 0
        file_zero_valid_events = 0
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            event_id = parts[0]
            sample_vals = parts[1:]
            actual_cols = len(sample_vals)
            file_total_cells += actual_cols
            events_in_file += 1
            psi_vals = []
            for val in sample_vals:
                try:
                    v = float(val)
                    if not np.isnan(v):
                        psi_vals.append(v)
                except ValueError:
                    continue
            valid_count = len(psi_vals)
            file_event_valid_counts.append(valid_count)
            valid_cells_in_file += valid_count
            if valid_count > 0:
                event_mean = sum(psi_vals) / valid_count
                file_event_means.append(event_mean)
                file_psi_values.extend(psi_vals)
                valid_events_in_file += 1
            else:
                file_zero_valid_events += 1
        total_events_all += events_in_file
        total_cells_all += file_total_cells
        valid_cells_all += valid_cells_in_file
        zero_valid_events += file_zero_valid_events
        all_psi_values.extend(file_psi_values)
        event_mean_psis.extend(file_event_means)
        event_valid_counts.extend(file_event_valid_counts)
        valid_psi_ratio_file = valid_cells_in_file / file_total_cells if file_total_cells > 0 else 0
        type_stats = {
            "type": etype,
            "total_events": events_in_file,
            "valid_events": valid_events_in_file,
            "valid_psi_ratio": valid_psi_ratio_file,
            "mean_event_psi": np.mean(file_event_means) if file_event_means else None,
            "global_mean_psi": np.mean(file_psi_values) if file_psi_values else None
        }
        type_stats_list.append(type_stats)

    if total_events_all == 0:
        return None

    global_mean_psi = np.mean(all_psi_values) if all_psi_values else 0
    valid_psi_ratio = valid_cells_all / total_cells_all if total_cells_all > 0 else 0
    dist = [0, 0, 0, 0]
    for mean in event_mean_psis:
        if mean <= 0.25:
            dist[0] += 1
        elif mean <= 0.5:
            dist[1] += 1
        elif mean <= 0.75:
            dist[2] += 1
        else:
            dist[3] += 1

    valid_counts_dist = defaultdict(int)
    for vc in event_valid_counts:
        valid_counts_dist[vc] += 1

    thresholds = [1, 2, 4, 6, 8]
    threshold_counts = {}
    for thresh in thresholds:
        count = sum(1 for vc in event_valid_counts if vc >= thresh)
        threshold_counts[str(thresh)] = count

    result = {
        "global_mean_psi": global_mean_psi,
        "total_psi_events": total_events_all,
        "total_psi_files": len(psi_files),
        "avg_samples_per_event": total_cells_all / total_events_all if total_events_all else 0,
        "valid_psi_count": int(valid_cells_all),
        "total_psi_cells": int(total_cells_all),
        "valid_psi_ratio": valid_psi_ratio,
        "psi_distribution": dist,
        "event_valid_sample_distribution": dict(valid_counts_dist),
        "threshold_counts": threshold_counts,
        "type_stats": type_stats_list,
        "zero_valid_events": zero_valid_events,
    }
    return result

# ------------------------------------------------------------
# Parse .ioe files and extract detailed event metadata
def parse_ioe_detailed(ioe_dir):
    """
    SUPPA .ioe file format (tab-separated):
    First line header, then each line: event_id [chrom start end strand] [extra columns]
    But the event_id may contain coordinates (e.g., "ENSG...;SE:chr1:1082987-1083915:1084086-1084353:-").
    We will try to extract coordinates from the event_id if present.
    """
    event_info = {}
    ioe_files = glob.glob(os.path.join(ioe_dir, "*.ioe"))
    for ioe_file in ioe_files:
        basename = os.path.basename(ioe_file)
        # Extract event type from filename
        match = re.search(r'[._]?events?_([A-Z0-9_]+)\.ioe', basename)
        if match:
            etype = match.group(1)
        else:
            etype = os.path.splitext(basename)[0]
            if etype.startswith('events_'):
                etype = etype[7:]
            if not etype:
                etype = "unknown"

        with open(ioe_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            try:
                header = next(reader)   # skip header
            except StopIteration:
                continue
            for row in reader:
                if not row:
                    continue
                event_id = row[0]

                # Attempt to parse coordinates from event_id if it contains ';'
                chrom = start = end = strand = ''
                extra = ''
                if ';' in event_id:
                    # Example: "ENSG00000131591.18;SE:chr1:1082987-1083915:1084086-1084353:-"
                    parts = event_id.split(';', 1)
                    # parts[0] is the gene/transcript part, parts[1] is event spec
                    event_id_clean = parts[0]
                    event_spec = parts[1] if len(parts) > 1 else ''
                    # Parse event_spec: e.g., "SE:chr1:1082987-1083915:1084086-1084353:-"
                    # For different event types, the format may differ
                    # We'll split by colon and try to extract chromosome, start, end, strand
                    spec_parts = event_spec.split(':')
                    if len(spec_parts) >= 5:
                        # Format: type:chrom:start-end:start-end:strand? Actually for SE it's type:chrom:start1-end1:start2-end2:strand
                        # We'll just take the first start and last end
                        try:
                            chrom = spec_parts[1]
                            # First interval
                            interval1 = spec_parts[2].split('-')
                            interval2 = spec_parts[3].split('-')
                            if len(interval1) == 2 and len(interval2) == 2:
                                start = min(int(interval1[0]), int(interval2[0]))
                                end = max(int(interval1[1]), int(interval2[1]))
                            else:
                                start = int(interval1[0])
                                end = int(interval1[1])
                            strand = spec_parts[-1] if len(spec_parts) > 4 else ''
                        except (ValueError, IndexError):
                            pass
                    extra = event_spec
                else:
                    # If no semicolon, fall back to using row columns if present
                    if len(row) >= 5:
                        try:
                            chrom = row[1]
                            start = int(row[2])
                            end = int(row[3])
                            strand = row[4]
                        except (ValueError, IndexError):
                            pass
                    extra = '\t'.join(row[5:]) if len(row) > 5 else ''

                # Store metadata
                event_info[event_id] = {
                    'type': etype,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'extra': extra
                }
    return event_info

# ------------------------------------------------------------
# Parse enhanced_transcript_tpm_all.tsv (may not exist)
def parse_tpm_file(tpm_path):
    if not os.path.exists(tpm_path):
        return None
    transcript_means = []
    with open(tpm_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            if len(row) < 9:
                continue
            tpm_values = [float(x) for x in row[1:9] if x]
            if tpm_values:
                transcript_means.append(sum(tpm_values) / len(tpm_values))
    total = len(transcript_means)
    avg_tpm = sum(transcript_means) / total if total else 0
    return {"total_transcripts": total, "mean_tpm": avg_tpm}

def parse_novel_file(novel_path):
    if not os.path.exists(novel_path):
        return None
    avg_tpms = []
    with open(novel_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                avg_tpms.append(float(row['avg_TPM']))
            except (KeyError, ValueError):
                continue
    total = len(avg_tpms)
    mean_avg_tpm = sum(avg_tpms) / total if total else 0
    return {"total_novel": total, "mean_avg_tpm": mean_avg_tpm}

# ------------------------------------------------------------
# Process a single directory and return statistics dict
def process_group(group_dir, group_name):
    print(f"Processing {group_name} group: {group_dir}")
    if not os.path.isdir(group_dir):
        raise ValueError(f"Directory not found: {group_dir}")

    events_dir = os.path.join(group_dir, "events")
    psi_dir = os.path.join(group_dir, "psi")
    tpm_file = os.path.join(group_dir, "enhanced_transcript_tpm_all.tsv")
    novel_file = os.path.join(group_dir, "high_quality_novel_transcripts.tsv")

    stats = {
        "group": group_name,
        "total_events": 0,
        "event_types": {},
        "psi_stats": None,
        "transcript_stats": None,
        "novel_stats": None,
        "psi_dir": psi_dir,
        "group_dir": group_dir
    }

    if os.path.isdir(events_dir):
        total, etypes = parse_ioe_files(events_dir)
        stats["total_events"] = total
        stats["event_types"] = etypes
    else:
        print(f"Warning: events folder not found in {group_dir}")

    if os.path.isdir(psi_dir):
        stats["psi_stats"] = parse_psi_files(psi_dir)
    else:
        print(f"Warning: psi folder not found in {group_dir}")

    stats["transcript_stats"] = parse_tpm_file(tpm_file)
    stats["novel_stats"] = parse_novel_file(novel_file)

    return stats

# ------------------------------------------------------------
# Save detailed event PSI data to SQLite
def save_detailed_psi_to_db(stats, tpm_matrix_path, db_path="stats.db"):
    """
    将 enhanced 组的详细 PSI 数据保存到 SQLite 数据库。
    只创建 events 和 event_psi 两张表，不涉及 groups 和 samples。
    events 表包含事件 ID、事件类型和所属组（enhanced）。
    event_psi 表存储每个事件在每个样本中的 PSI 值。
    """
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # 删除旧表（如果有）
    c.execute("DROP TABLE IF EXISTS events")
    c.execute("DROP TABLE IF EXISTS event_psi")

    # 创建 events 表
    c.execute('''CREATE TABLE events (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    event_id TEXT NOT NULL,
                    event_type TEXT NOT NULL,
                    group_name TEXT
                )''')

    # 创建 event_psi 表
    c.execute('''CREATE TABLE event_psi (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    event_id INTEGER,
                    sample_name TEXT,
                    psi_value REAL,
                    FOREIGN KEY(event_id) REFERENCES events(id)
                )''')

    # 读取样本名（从 TPM 矩阵的第一行）
    if not os.path.exists(tpm_matrix_path):
        print(f"Warning: TPM matrix not found at {tpm_matrix_path}, using generic sample names")
        sample_names = None
    else:
        with open(tpm_matrix_path, 'r') as f:
            first_line = f.readline().strip()
            sample_names = first_line.split('\t')
        print(f"Loaded {len(sample_names)} sample names from TPM matrix")

    # 确定样本名（如果没有，从 PSI 文件推测）
    psi_dir = stats.get("psi_dir")
    if not psi_dir or not os.path.isdir(psi_dir):
        print(f"Warning: PSI directory not found: {psi_dir}")
        conn.close()
        return

    psi_files = glob.glob(os.path.join(psi_dir, "*.psi"))
    if not psi_files:
        print(f"Warning: No PSI files found in {psi_dir}")
        conn.close()
        return

    # 如果 sample_names 为空，从第一个 PSI 文件推测
    if not sample_names:
        with open(psi_files[0], 'r') as f:
            first_line = f.readline().strip()
            parts = first_line.split('\t')
            num_samples = len(parts) - 1
            sample_names = [f"sample_{i+1}" for i in range(num_samples)]
        print(f"Using generic sample names: {num_samples} samples")

    # 遍历所有 PSI 文件
    for psi_file in psi_files:
        basename = os.path.basename(psi_file)
        # 从文件名提取事件类型，例如 "SE_SE_strict.psi" -> "SE"
        match = re.match(r'([A-Z0-9]+)(?:_[A-Z0-9]+)?\.psi', basename)
        if match:
            event_type_from_file = match.group(1)
        else:
            event_type_from_file = "unknown"

        with open(psi_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                event_id = parts[0]
                sample_vals = parts[1:]

                # 确定最终的事件类型：优先从事件 ID 中提取（如果有 "SE:" 等），否则用文件名提取的类型
                event_type = event_type_from_file
                if ':' in event_id:
                    # 尝试从事件 ID 中提取类型，格式如 "ENSG...;SE:..." 或 "MSTRG.123;SE:..."
                    id_parts = event_id.split(';')
                    if len(id_parts) > 1:
                        spec = id_parts[1]
                        if ':' in spec:
                            event_type = spec.split(':')[0]   # 取冒号前的部分，如 "SE"
                # 如果还是 unknown，保留文件名提取的类型

                # 插入 events 表（如果已存在，不重复插入，但这里简单起见，每次都尝试插入，利用唯一性约束）
                c.execute("INSERT OR IGNORE INTO events (event_id, event_type, group_name) VALUES (?,?,?)",
                          (event_id, event_type, "enhanced"))
                c.execute("SELECT id FROM events WHERE event_id = ?", (event_id,))
                event_db_id = c.fetchone()[0]

                # 插入每个样本的 PSI 值
                for sample_idx, val in enumerate(sample_vals):
                    if sample_idx >= len(sample_names):
                        continue
                    try:
                        psi_val = float(val)
                        if not np.isnan(psi_val):
                            c.execute("INSERT INTO event_psi (event_id, sample_name, psi_value) VALUES (?,?,?)",
                                      (event_db_id, sample_names[sample_idx], psi_val))
                    except ValueError:
                        continue

    conn.commit()
    conn.close()
    print(f"Detailed PSI data (enhanced group) saved to {db_path}")

# ------------------------------------------------------------
# Generate HTML report (unchanged, kept as before)
def generate_html(stats_list, output="index.html"):
    stats_js = json.dumps(stats_list, indent=2, default=str)
    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AS Event Statistics Comparison · Pie Charts</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <style>
        body {{
            background-color: #0d1117;
            color: #c9d1d9;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
            line-height: 1.5;
            padding: 2rem;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        h1 {{
            border-bottom: 1px solid #30363d;
            padding-bottom: 0.5rem;
        }}
        .selector {{
            margin: 1rem 0;
        }}
        select {{
            background-color: #21262d;
            color: #c9d1d9;
            border: 1px solid #30363d;
            padding: 0.5rem 1rem;
            border-radius: 6px;
            font-size: 1rem;
        }}
        .stats-card {{
            background-color: #161b22;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 1.5rem;
            margin-top: 1rem;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 1rem;
        }}
        .stat-item {{
            background-color: #21262d;
            padding: 1rem;
            border-radius: 6px;
        }}
        .stat-label {{
            font-size: 0.9rem;
            color: #8b949e;
        }}
        .stat-value {{
            font-size: 2rem;
            font-weight: 600;
        }}
        .chart-wrapper {{
            display: flex;
            gap: 20px;
            margin-top: 20px;
            flex-wrap: nowrap;
        }}
        .chart-container {{
            flex: 2 1 0;
            min-width: 0;
            background-color: #0d1117;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 20px;
            height: 300px;
        }}
        .legend-panel {{
            flex: 1 1 250px;
            background-color: #0d1117;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 16px;
            overflow-y: auto;
            max-height: 300px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin-bottom: 8px;
            padding: 4px 8px;
            border-radius: 4px;
            transition: background-color 0.2s;
        }}
        .legend-item:hover {{
            background-color: #21262d;
        }}
        .color-box {{
            display: inline-block;
            width: 16px;
            height: 16px;
            border-radius: 4px;
            margin-right: 8px;
            flex-shrink: 0;
        }}
        .legend-desc {{
            flex: 1;
            font-size: 13px;
            margin-right: 8px;
            word-break: break-word;
            color: #c9d1d9;
        }}
        .legend-count {{
            color: #8b949e;
            font-family: monospace;
            flex-shrink: 0;
        }}
        h3, h4 {{
            color: #c9d1d9;
            border-bottom: 1px solid #30363d;
            padding-bottom: 0.3rem;
            margin-top: 1.5rem;
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>AS Event Statistics Comparison</h1>
    <div class="selector">
        <label for="group-select">Select annotation group: </label>
        <select id="group-select">
            <option value="expressed" selected>expressed (gencode.v45.expressed.gtf)</option>
            <option value="enhanced">enhanced (gencode.v45.enhanced.gtf)</option>
        </select>
    </div>

    <div id="stats-display">
        <!-- dynamic content -->
    </div>
</div>

<script>
    const statsData = {stats_js};

    let pieTypeChart = null;
    let pieSampleChart = null;
    const colorPalette = [
        '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0', '#9966FF', '#FF9F40',
        '#E7E9ED', '#76A346', '#F7464A', '#46BFBD', '#FDB45C', '#949FB1',
        '#4D5360', '#1FC8F8', '#9AC0CD', '#FF8A80', '#B39DDB', '#80CBC4'
    ];

    function buildLegendHTML(labels, counts, colors) {{
        let legend = '';
        for (let i = 0; i < labels.length; i++) {{
            legend += `
                <div class="legend-item">
                    <span class="color-box" style="background-color: ${{colors[i]}};"></span>
                    <span class="legend-desc">${{labels[i]}}</span>
                    <span class="legend-count">${{counts[i].toLocaleString()}}</span>
                </div>
            `;
        }}
        return legend;
    }}

    function renderGroup(groupName) {{
        const group = statsData.find(g => g.group === groupName);
        if (!group) return '<div class="stats-card">No data available</div>';

        let html = `<div class="stats-card">`;
        html += `<div class="stat-item"><span class="stat-label">Total events</span><div class="stat-value">${{group.total_events}}</div></div>`;
        html += `<h3>Event type distribution</h3>`;
        html += `<div class="chart-wrapper">`;
        html += `<div class="chart-container"><canvas id="pie-event-types"></canvas></div>`;
        html += `<div class="legend-panel" id="legend-event-types"></div>`;
        html += `</div>`;

        if (group.psi_stats) {{
            html += `<h3>PSI statistics</h3>`;
            html += `<div class="stats-grid">`;
            html += `<div class="stat-item"><span class="stat-label">Global mean PSI</span><div class="stat-value">${{group.psi_stats.global_mean_psi.toFixed(3)}}</div></div>`;
            html += `<div class="stat-item"><span class="stat-label">Events with PSI</span><div class="stat-value">${{group.psi_stats.total_psi_events}}</div></div>`;
            html += `<div class="stat-item"><span class="stat-label">Valid PSI ratio</span><div class="stat-value">${{(group.psi_stats.valid_psi_ratio * 100).toFixed(1)}}%</div></div>`;
            html += `</div>`;
            html += `<h4>Events by valid sample count</h4>`;
            html += `<div class="chart-wrapper">`;
            html += `<div class="chart-container"><canvas id="pie-sample-dist"></canvas></div>`;
            html += `<div class="legend-panel" id="legend-sample-dist"></div>`;
            html += `</div>`;
        }}

        if (group.transcript_stats) {{
            html += `<h3>All transcripts (TPM)</h3>`;
            html += `<div class="stats-grid">`;
            html += `<div class="stat-item"><span class="stat-label">Total transcripts</span><div class="stat-value">${{group.transcript_stats.total_transcripts}}</div></div>`;
            html += `<div class="stat-item"><span class="stat-label">Mean TPM</span><div class="stat-value">${{group.transcript_stats.mean_tpm.toFixed(2)}}</div></div>`;
            html += `</div>`;
        }}

        if (group.novel_stats) {{
            html += `<h3>High-quality novel transcripts</h3>`;
            html += `<div class="stats-grid">`;
            html += `<div class="stat-item"><span class="stat-label">Novel transcripts</span><div class="stat-value">${{group.novel_stats.total_novel}}</div></div>`;
            html += `<div class="stat-item"><span class="stat-label">Mean avg_TPM</span><div class="stat-value">${{group.novel_stats.mean_avg_tpm.toFixed(2)}}</div></div>`;
            html += `</div>`;
        }}

        html += `</div>`;
        return html;
    }}

    function drawCharts(groupName) {{
        const group = statsData.find(g => g.group === groupName);
        if (!group) return;

        if (pieTypeChart) pieTypeChart.destroy();
        if (pieSampleChart) pieSampleChart.destroy();

        const eventTypes = group.event_types || {{}};
        const typeLabels = Object.keys(eventTypes);
        const typeCounts = Object.values(eventTypes);
        const typeColors = typeLabels.map((_, i) => colorPalette[i % colorPalette.length]);

        const legendTypeDiv = document.getElementById('legend-event-types');
        if (legendTypeDiv) {{
            legendTypeDiv.innerHTML = buildLegendHTML(typeLabels, typeCounts, typeColors);
        }}

        const ctxType = document.getElementById('pie-event-types')?.getContext('2d');
        if (ctxType) {{
            pieTypeChart = new Chart(ctxType, {{
                type: 'pie',
                data: {{
                    labels: typeLabels,
                    datasets: [{{
                        data: typeCounts,
                        backgroundColor: typeColors,
                        borderColor: '#0d1117',
                        borderWidth: 2
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{ display: false }},
                        tooltip: {{ callbacks: {{ label: (ctx) => `${{ctx.label}}: ${{ctx.raw.toLocaleString()}}` }} }}
                    }}
                }}
            }});
        }}

        const psiStats = group.psi_stats;
        if (psiStats) {{
            const dist = psiStats.event_valid_sample_distribution || {{}};
            const totalEvents = group.total_events;
            const totalPSI = psiStats.total_psi_events;
            const zeroValid = psiStats.zero_valid_events || 0;
            const noPSI = totalEvents - totalPSI;

            const sampleKeys = Object.keys(dist).map(Number);
            const maxSample = Math.max(...sampleKeys, 0);

            const sampleLabels = [];
            const sampleCounts = [];

            for (let i = 0; i <= maxSample; i++) {{
                const cnt = dist[i] || 0;
                if (cnt > 0) {{
                    sampleLabels.push(`${{i}} sample${{i !== 1 ? 's' : ''}}`);
                    sampleCounts.push(cnt);
                }}
            }}

            if (zeroValid > 0) {{
                sampleLabels.push('0 samples (all NaN)');
                sampleCounts.push(zeroValid);
            }}

            if (noPSI > 0) {{
                sampleLabels.push('no PSI data');
                sampleCounts.push(noPSI);
            }}

            const sampleColors = sampleLabels.map((_, idx) => colorPalette[(idx + 5) % colorPalette.length]);

            const legendSampleDiv = document.getElementById('legend-sample-dist');
            if (legendSampleDiv) {{
                legendSampleDiv.innerHTML = buildLegendHTML(sampleLabels, sampleCounts, sampleColors);
            }}

            const ctxSample = document.getElementById('pie-sample-dist')?.getContext('2d');
            if (ctxSample) {{
                pieSampleChart = new Chart(ctxSample, {{
                    type: 'pie',
                    data: {{
                        labels: sampleLabels,
                        datasets: [{{
                            data: sampleCounts,
                            backgroundColor: sampleColors,
                            borderColor: '#0d1117',
                            borderWidth: 2
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {{
                            legend: {{ display: false }},
                            tooltip: {{ callbacks: {{ label: (ctx) => `${{ctx.label}}: ${{ctx.raw.toLocaleString()}}` }} }}
                        }}
                    }}
                }});
            }}
        }}
    }}

    function updateDisplay() {{
        const select = document.getElementById('group-select');
        const groupName = select.value;
        const displayDiv = document.getElementById('stats-display');
        displayDiv.innerHTML = renderGroup(groupName);
        drawCharts(groupName);
    }}

    window.onload = function() {{
        updateDisplay();
        document.getElementById('group-select').addEventListener('change', updateDisplay);
    }};
</script>
</body>
</html>
"""
    with open(output, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print(f"HTML generated: {output}")

def main():
    parser = argparse.ArgumentParser(description="Process AS data and generate statistics report")
    parser.add_argument("dirs", nargs=2, help="Two directory paths, in order: expressed, enhanced")
    args = parser.parse_args()

    group_names = ["expressed", "enhanced"]
    stats_list = []
    for i, dir_path in enumerate(args.dirs):
        try:
            stats = process_group(dir_path, group_names[i])
            stats_list.append(stats)
        except Exception as e:
            print(f"Error processing {dir_path}: {e}")
            sys.exit(1)

    # For standalone usage, need TPM matrix path; but in pipeline we call with extra arg.
    # We'll keep as is and expect it to be passed via pipeline call.
    # In pipeline, save_detailed_psi_to_db is called with tpm_matrix_path argument.
    # Here for direct call, we skip DB creation or use default.
    # We'll just generate HTML.
    generate_html(stats_list, "index.html")
    print("Done.")

if __name__ == "__main__":
    main()