#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
APA 分析完整流水线（全前端图表版）
用法: python run_apa_pipeline.py --all_results all_apa_results.txt --predictions apa_predictions.txt \
        --utr refined_3utr.bed --genome hg38.fa --outdir ./output --prefix my_analysis
"""

import os
import argparse
import sys
import subprocess
import sqlite3
import json
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist


# ---------- 1. apa_summary 模块 ----------
def apa_summary(all_results_file, outdir, prefix=None):
    """生成位点级别汇总表"""
    print("正在执行 apa_summary 模块...")
    df = pd.read_csv(all_results_file, sep='\t')

    sep_cols = [col for col in df.columns if 'Separate_Exp' in col]
    total_cols = [col for col in df.columns if 'Total_Exp' in col]

    total_dict = {}
    for idx, row in df.iterrows():
        gene = row['Gene']
        for col in total_cols:
            sample = col.replace('Group_1_', '').replace('_Total_Exp', '')
            total_dict[(gene, sample)] = row[col]

    records = []
    for idx, row in df.iterrows():
        gene = row['Gene']
        loci = row['Loci']
        proximal_sites = [] if pd.isna(row['Predicted_APA']) or row['Predicted_APA'] == '' else row[
            'Predicted_APA'].split(',')
        n_total_sites = len(proximal_sites) + 1

        try:
            distal_coord = loci.split(':')[1].split('-')[1]
        except IndexError:
            print(f"警告: 基因 {gene} 的 Loci 格式异常: {loci}")
            distal_coord = 'unknown'

        for col in sep_cols:
            sample = col.replace('Group_1_', '').replace('_Separate_Exp', '')
            exp_str = row[col]
            if pd.isna(exp_str) or exp_str == '':
                continue
            exp_values = list(map(float, exp_str.split(',')))
            if len(exp_values) != n_total_sites:
                continue

            for site_idx, exp in enumerate(exp_values):
                coord = proximal_sites[site_idx] if site_idx < len(proximal_sites) else distal_coord
                records.append({
                    'Gene': gene,
                    'Loci': loci,
                    'Site_Index': site_idx + 1,
                    'Site_Coord': coord,
                    'Sample': sample,
                    'Expression': exp,
                    'Total_Exp': total_dict.get((gene, sample), np.nan)
                })

    long_df = pd.DataFrame(records)
    if long_df.empty:
        print("错误：没有有效数据，请检查输入文件格式。")
        sys.exit(1)

    long_df['Rel_Usage'] = long_df['Expression'] / long_df['Total_Exp']

    summary = long_df.groupby(['Gene', 'Loci', 'Site_Index', 'Site_Coord']).agg(
        mean_expression=('Expression', 'mean'),
        std_expression=('Expression', 'std'),
        mean_rel_usage=('Rel_Usage', 'mean'),
        std_rel_usage=('Rel_Usage', 'std'),
        n_samples=('Sample', 'count')
    ).reset_index()

    out_file = os.path.join(outdir, f"{prefix}_apa_summary.txt" if prefix else "apa_summary.txt")
    summary.to_csv(out_file, sep='\t', index=False)
    print(f"apa_summary 完成: {out_file}")
    return summary, long_df


# ---------- 2. deapa 模块 ----------
def generate_apa_sites(pred_df, output_bed):
    records = []
    for idx, row in pred_df.iterrows():
        chrom = row['chrom']
        sites = str(row['Predicted_APA']).split(',')
        for i, site in enumerate(sites):
            site = site.strip()
            coord = site.split(':')[1].split('-')[0] if ':' in site else site
            try:
                start = int(coord) - 1
                end = start + 1
            except:
                continue
            records.append([chrom, start, end, f"{row['Gene']}_site{i + 1}"])
    bed_df = pd.DataFrame(records, columns=['chrom', 'start', 'end', 'name'])
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"生成 {output_bed}，共 {len(bed_df)} 个位点")


def generate_apa_sequences(apa_sites_bed, genome_fa, output_fa):
    cmd = f"bedtools getfasta -fi {genome_fa} -bed {apa_sites_bed} -fo {output_fa}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"生成 {output_fa}")


def generate_gene_level_summary(pred_df, output_file):
    summary = []
    for gene, group in pred_df.groupby('Gene'):
        all_sites = set()
        for sites in group['Predicted_APA'].dropna():
            for s in str(sites).split(','):
                all_sites.add(s.strip())
        num_sites = len(all_sites)
        summary.append({'Gene': gene, 'num_apa_sites': num_sites, 'has_multiple_apa': num_sites > 1})
    pd.DataFrame(summary).to_csv(output_file, sep='\t', index=False)
    print(f"生成 {output_file}")


def generate_detailed_analysis(pred_df, output_file):
    exp_cols = [col for col in pred_df.columns if 'Separate_Exp' in col]
    if not exp_cols:
        detailed = pred_df[['Gene', 'chrom', 'Predicted_APA']].copy()
        detailed['Mean_Usage'] = np.nan
    else:
        records = []
        for idx, row in pred_df.iterrows():
            chrom = row['chrom']
            sites = [s.strip() for s in str(row['Predicted_APA']).split(',')]
            usages = defaultdict(list)
            for col in exp_cols:
                vals = str(row[col]).split(',')
                for i, val in enumerate(vals):
                    if i < len(sites):
                        try:
                            usages[sites[i]].append(float(val))
                        except:
                            pass
            for site, usage_list in usages.items():
                if usage_list:
                    records.append(
                        {'Gene': row['Gene'], 'chrom': chrom, 'APA_site': site, 'Mean_Usage': np.mean(usage_list)})
        detailed = pd.DataFrame(records)
    detailed.to_csv(output_file, sep='\t', index=False)
    print(f"生成 {output_file}")


def generate_novel_transcripts_summary(utr_bed, output_file):
    if not os.path.exists(utr_bed):
        pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'score', 'strand']).to_csv(output_file, sep='\t',
                                                                                          index=False)
        return
    df = pd.read_csv(utr_bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    df.to_csv(output_file, sep='\t', index=False)
    print(f"生成 {output_file}")


def deapa_module(predictions_file, utr_file, genome_file, outdir, prefix=None):
    print("正在执行 deapa 模块...")
    pred_df = pd.read_csv(predictions_file, sep='\t')
    required_cols = ['Gene', 'Predicted_APA', 'Loci']
    for col in required_cols:
        if col not in pred_df.columns:
            raise ValueError(f"缺少必要列: {col}")

    def extract_chrom(gene):
        parts = str(gene).split('|')
        return parts[2] if len(parts) >= 3 else 'unknown'

    pred_df['chrom'] = pred_df['Gene'].apply(extract_chrom)

    sites_bed = os.path.join(outdir, f"{prefix}_apa_sites.bed" if prefix else "apa_sites.bed")
    generate_apa_sites(pred_df, sites_bed)

    seq_fa = os.path.join(outdir, f"{prefix}_apa_sequences.fa" if prefix else "apa_sequences.fa")
    generate_apa_sequences(sites_bed, genome_file, seq_fa)

    gene_summary = os.path.join(outdir,
                                f"{prefix}_apa_gene_level_summary.txt" if prefix else "apa_gene_level_summary.txt")
    generate_gene_level_summary(pred_df, gene_summary)

    detailed = os.path.join(outdir, f"{prefix}_apa_detailed_analysis.txt" if prefix else "apa_detailed_analysis.txt")
    generate_detailed_analysis(pred_df, detailed)

    novel = os.path.join(outdir,
                         f"{prefix}_apa_novel_transcripts_summary.txt" if prefix else "apa_novel_transcripts_summary.txt")
    generate_novel_transcripts_summary(utr_file, novel)

    print("deapa 模块完成。")


# ---------- 3. further 模块（仅计算数据，不绘图）----------
def further_module(all_results_file, outdir, prefix):
    """计算所有图表所需的数据并保存为 JSON"""
    print("正在执行 further 模块（计算前端绘图数据）...")
    df = pd.read_csv(all_results_file, sep='\t')

    sep_cols = [col for col in df.columns if 'Separate_Exp' in col]
    total_cols = [col for col in df.columns if 'Total_Exp' in col]

    total_dict = {}
    for idx, row in df.iterrows():
        gene = row['Gene']
        for col in total_cols:
            sample = col.replace('Group_1_', '').replace('_Total_Exp', '')
            total_dict[(gene, sample)] = row[col]

    records = []
    for idx, row in df.iterrows():
        gene = row['Gene']
        loci = row['Loci']
        proximal_sites = [] if pd.isna(row['Predicted_APA']) or row['Predicted_APA'] == '' else row[
            'Predicted_APA'].split(',')
        n_total_sites = len(proximal_sites) + 1
        try:
            distal_coord = loci.split(':')[1].split('-')[1]
        except IndexError:
            distal_coord = 'unknown'

        for col in sep_cols:
            sample = col.replace('Group_1_', '').replace('_Separate_Exp', '')
            exp_str = row[col]
            if pd.isna(exp_str) or exp_str == '':
                continue
            exp_values = list(map(float, exp_str.split(',')))
            if len(exp_values) != n_total_sites:
                continue
            for site_idx, exp in enumerate(exp_values):
                coord = proximal_sites[site_idx] if site_idx < len(proximal_sites) else distal_coord
                records.append({
                    'Gene': gene, 'Loci': loci, 'Site_Index': site_idx + 1, 'Site_Coord': coord,
                    'Sample': sample, 'Expression': exp, 'Total_Exp': total_dict.get((gene, sample), np.nan)
                })

    long_df = pd.DataFrame(records)
    long_df['Rel_Usage'] = long_df['Expression'] / long_df['Total_Exp']

    summary = long_df.groupby(['Gene', 'Loci', 'Site_Index', 'Site_Coord']).agg(
        mean_expression=('Expression', 'mean'),
        std_expression=('Expression', 'std'),
        mean_rel_usage=('Rel_Usage', 'mean'),
        std_rel_usage=('Rel_Usage', 'std'),
        n_samples=('Sample', 'count')
    ).reset_index()

    # 保存汇总表
    summary.to_csv(os.path.join(outdir, f"{prefix}_summary.txt"), sep='\t', index=False)

    # ----- 准备前端绘图数据 -----
    plot_data = {}

    # 1. 位点数量分布
    site_count = summary.groupby('Gene')['Site_Index'].max().value_counts().sort_index()
    plot_data['site_count_dist'] = {
        'labels': [int(x) for x in site_count.index],
        'values': site_count.values.tolist()
    }

    # 2. 主要位点使用率分布
    major_site = summary.loc[summary.groupby('Gene')['mean_rel_usage'].idxmax()]
    plot_data['major_usage_dist'] = major_site['mean_rel_usage'].tolist()

    # 3. CV 数据
    df_cv = summary[summary['mean_expression'] > 0].copy()
    df_cv['cv_expression'] = df_cv['std_expression'] / df_cv['mean_expression']
    df_cv['cv_rel_usage'] = df_cv['std_rel_usage'] / df_cv['mean_rel_usage']
    plot_data['cv_expression'] = df_cv['cv_expression'].dropna().tolist()
    plot_data['cv_rel_usage'] = df_cv['cv_rel_usage'].dropna().tolist()

    # 4. 宽表（用于热图与聚类）
    wide_df = summary.pivot_table(index='Gene', columns='Site_Index', values='mean_rel_usage', aggfunc='first').fillna(
        0)
    wide_df.columns = [f'Site_{int(col)}' for col in wide_df.columns]
    wide_df = wide_df.loc[(wide_df > 0).sum(axis=1) >= 2]  # 至少两个位点

    if len(wide_df) > 1:
        # 层次聚类
        Z = linkage(wide_df.values, method='ward', metric='euclidean')
        order = leaves_list(Z)
        ordered_genes = wide_df.index[order].tolist()
        plot_data['dendrogram_order'] = ordered_genes
        plot_data['heatmap_data'] = {
            'genes': ordered_genes,
            'columns': wide_df.columns.tolist(),
            'values': wide_df.loc[ordered_genes].values.tolist()
        }

        # t-SNE 降维坐标（若 sklearn 可用）
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            data_scaled = StandardScaler().fit_transform(wide_df.values)
            perplexity = min(50, len(wide_df) - 1)
            tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
            coords = tsne.fit_transform(data_scaled)
            clusters = fcluster(Z, t=10, criterion='maxclust')
            plot_data['tsne'] = {
                'x': coords[:, 0].tolist(),
                'y': coords[:, 1].tolist(),
                'genes': wide_df.index.tolist(),
                'clusters': clusters.tolist()
            }
        except ImportError:
            plot_data['tsne'] = None
    else:
        plot_data['dendrogram_order'] = []
        plot_data['heatmap_data'] = None
        plot_data['tsne'] = None

    # 保存 JSON
    json_path = os.path.join(outdir, f"{prefix}_plot_data.json")
    with open(json_path, 'w') as f:
        json.dump(plot_data, f, indent=2)
    print(f"前端绘图数据已保存: {json_path}")


# ---------- 4. HTML 报告生成（全前端渲染）----------
def extract_chromosome(gene_str):
    if pd.isna(gene_str):
        return 'unknown'
    parts = str(gene_str).split('|')
    return parts[2] if len(parts) >= 3 else 'unknown'


def generate_html_report(outdir, prefix, pred_file, utr_file, plot_data_json):
    """生成包含 Plotly 交互图表的深色风格报告"""
    # 读取必要文件
    gene_summary = pd.read_csv(os.path.join(outdir, f"{prefix}_apa_gene_level_summary.txt"), sep='\t')
    detailed = pd.read_csv(os.path.join(outdir, f"{prefix}_apa_detailed_analysis.txt"), sep='\t')
    sites_bed = os.path.join(outdir, f"{prefix}_apa_sites.bed")
    apa_sites = pd.read_csv(sites_bed, sep='\t', header=None,
                            names=['chrom', 'start', 'end', 'name']) if os.path.exists(sites_bed) else pd.DataFrame()

    # 染色体分布数据
    gene_chrom_counts = gene_summary['Gene'].apply(extract_chromosome).value_counts().to_dict()
    site_chrom_counts = detailed['Gene'].apply(extract_chromosome).value_counts().to_dict()
    detailed['Mean_Usage'] = pd.to_numeric(detailed['Mean_Usage'], errors='coerce')
    chrom_usage = {}
    for chrom in site_chrom_counts.keys():
        mask = detailed['Gene'].apply(extract_chromosome) == chrom
        usages = detailed.loc[mask, 'Mean_Usage'].dropna().tolist()
        if usages:
            chrom_usage[chrom] = usages

    # 统计卡片数据
    stats = {
        'total_genes': len(gene_summary),
        'total_apa_sites': len(apa_sites),
        'multi_apa_genes': int(
            gene_summary['has_multiple_apa'].sum()) if 'has_multiple_apa' in gene_summary.columns else 0,
        'mean_usage_mean': detailed['Mean_Usage'].mean() if not detailed.empty else 0
    }

    # 读取前端绘图数据
    with open(plot_data_json, 'r') as f:
        plot_data = json.load(f)

    # 构建 HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>APA Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{
            background-color: #0d1117;
            color: #c9d1d9;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
            padding: 20px;
            margin: 0;
        }}
        h1, h2, h3 {{
            color: #f0f6fc;
            border-bottom: 1px solid #30363d;
            padding-bottom: 0.3em;
        }}
        .summary-cards {{
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            margin-bottom: 30px;
        }}
        .card {{
            background-color: #161b22;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 20px;
            min-width: 180px;
            flex: 1 0 auto;
        }}
        .card h3 {{
            margin-top: 0;
            border-bottom: none;
            font-size: 1.1em;
            color: #8b949e;
        }}
        .card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #79c0ff;
        }}
        .chart-container {{
            background-color: #161b22;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 20px;
            margin-bottom: 30px;
        }}
        select, button {{
            background-color: #0d1117;
            color: #c9d1d9;
            border: 1px solid #30363d;
            border-radius: 6px;
            padding: 8px 12px;
            margin-right: 10px;
        }}
        .tab {{
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
        }}
        .tab button {{
            background-color: #21262d;
            border: none;
        }}
        .tab button.active {{
            background-color: #30363d;
            color: #f0f6fc;
        }}
        .plotly-graph-div {{
            border-radius: 6px;
        }}
    </style>
</head>
<body>
    <h1>APA Analysis Report</h1>

    <div class="summary-cards">
        <div class="card"><h3>Total Genes</h3><div class="value">{stats['total_genes']}</div></div>
        <div class="card"><h3>Total APA Sites</h3><div class="value">{stats['total_apa_sites']}</div></div>
        <div class="card"><h3>Multi-APA Genes</h3><div class="value">{stats['multi_apa_genes']}</div></div>
        <div class="card"><h3>Mean Usage (avg)</h3><div class="value">{stats['mean_usage_mean']:.3f}</div></div>
    </div>

    <h2>Chromosome Distribution</h2>
    <div class="chart-container" id="chromosome-charts"></div>

    <h2>APA Site Usage</h2>
    <div class="control-panel">
        <label for="chrom-select">Chromosome: </label>
        <select id="chrom-select"></select>
    </div>
    <div class="chart-container" id="usage-histogram"></div>

    <h2>Advanced Analysis</h2>
    <div class="tab">
        <button class="tablinks active" onclick="showTab('site-count')">Site Count</button>
        <button class="tablinks" onclick="showTab('major-usage')">Major Usage</button>
        <button class="tablinks" onclick="showTab('cv')">CV</button>
        <button class="tablinks" onclick="showTab('heatmap')">Heatmap</button>
        <button class="tablinks" onclick="showTab('tsne')">t-SNE</button>
    </div>

    <div id="site-count" class="chart-container tab-content" style="display:block;"></div>
    <div id="major-usage" class="chart-container tab-content" style="display:none;"></div>
    <div id="cv" class="chart-container tab-content" style="display:none;"></div>
    <div id="heatmap" class="chart-container tab-content" style="display:none;"></div>
    <div id="tsne" class="chart-container tab-content" style="display:none;"></div>

    <script>
        // 深色主题模板
        const darkTemplate = {{
            layout: {{
                paper_bgcolor: '#161b22',
                plot_bgcolor: '#161b22',
                font: {{ color: '#c9d1d9' }},
                xaxis: {{ gridcolor: '#30363d', zerolinecolor: '#30363d' }},
                yaxis: {{ gridcolor: '#30363d', zerolinecolor: '#30363d' }}
            }}
        }};

        // 嵌入数据
        const geneChromCounts = {json.dumps(gene_chrom_counts)};
        const siteChromCounts = {json.dumps(site_chrom_counts)};
        const chromUsage = {json.dumps(chrom_usage)};
        const plotData = {json.dumps(plot_data)};

        // ---------- 染色体分布图 ----------
        function naturalSort(chroms) {{
            return chroms.sort((a,b) => {{
                const pa = a.match(/^chr(\\d+)$/i);
                const pb = b.match(/^chr(\\d+)$/i);
                if(pa && pb) return parseInt(pa[1]) - parseInt(pb[1]);
                if(pa) return -1;
                if(pb) return 1;
                return a.localeCompare(b);
            }});
        }}

        function drawChromosomeCharts() {{
            const geneChroms = naturalSort(Object.keys(geneChromCounts));
            const siteChroms = naturalSort(Object.keys(siteChromCounts));

            const trace1 = {{
                x: geneChroms,
                y: geneChroms.map(c => geneChromCounts[c]),
                type: 'bar',
                name: 'Genes',
                marker: {{ color: '#79c0ff' }}
            }};
            const trace2 = {{
                x: siteChroms,
                y: siteChroms.map(c => siteChromCounts[c]),
                type: 'bar',
                name: 'APA Sites',
                marker: {{ color: '#ff7b72' }}
            }};

            Plotly.newPlot('chromosome-charts', [trace1, trace2], {{
                ...darkTemplate.layout,
                title: 'Gene and APA Site Counts per Chromosome',
                barmode: 'group'
            }});
        }}

        // ---------- 使用率直方图（按染色体）----------
        function populateChromSelect() {{
            const select = document.getElementById('chrom-select');
            select.innerHTML = '<option value="all">All Chromosomes</option>';
            naturalSort(Object.keys(chromUsage)).forEach(chrom => {{
                const opt = document.createElement('option');
                opt.value = chrom;
                opt.textContent = chrom;
                select.appendChild(opt);
            }});
        }}

        function drawUsageHistogram(selectedChrom) {{
            let data = selectedChrom === 'all' 
                ? [].concat(...Object.values(chromUsage))
                : chromUsage[selectedChrom] || [];

            const trace = {{
                x: data,
                type: 'histogram',
                marker: {{ color: '#7ee07f' }},
                nbinsx: 20
            }};

            Plotly.newPlot('usage-histogram', [trace], {{
                ...darkTemplate.layout,
                title: `APA Site Usage Distribution (${{selectedChrom}})`,
                xaxis: {{ title: 'Mean Usage', gridcolor: '#30363d' }},
                yaxis: {{ title: 'Frequency', type: 'log', gridcolor: '#30363d' }}
            }});
        }}

        // ---------- 标签页切换 ----------
        function showTab(tabId) {{
            document.querySelectorAll('.tab-content').forEach(el => el.style.display = 'none');
            document.getElementById(tabId).style.display = 'block';
            document.querySelectorAll('.tablinks').forEach(btn => btn.classList.remove('active'));
            event.target.classList.add('active');
        }}

        // ---------- 绘制高级图表 ----------
        function drawSiteCountDist() {{
            const trace = {{
                x: plotData.site_count_dist.labels,
                y: plotData.site_count_dist.values,
                type: 'bar',
                marker: {{ color: '#a371f7' }}
            }};
            Plotly.newPlot('site-count', [trace], {{
                ...darkTemplate.layout,
                title: 'Number of APA Sites per Gene',
                xaxis: {{ title: 'Number of Sites', tickvals: plotData.site_count_dist.labels }},
                yaxis: {{ title: 'Frequency' }}
            }});
        }}

        function drawMajorUsageDist() {{
            const trace = {{
                x: plotData.major_usage_dist,
                type: 'histogram',
                marker: {{ color: '#f0883e' }},
                nbinsx: 30
            }};
            Plotly.newPlot('major-usage', [trace], {{
                ...darkTemplate.layout,
                title: 'Mean Relative Usage of Major APA Site',
                xaxis: {{ title: 'Mean Relative Usage' }},
                yaxis: {{ title: 'Number of Genes' }}
            }});
        }}

        function drawCVDist() {{
            const trace1 = {{
                x: plotData.cv_expression,
                type: 'histogram',
                name: 'Expression CV',
                marker: {{ color: '#79c0ff' }},
                opacity: 0.7,
                nbinsx: 30
            }};
            const trace2 = {{
                x: plotData.cv_rel_usage,
                type: 'histogram',
                name: 'Relative Usage CV',
                marker: {{ color: '#7ee07f' }},
                opacity: 0.7,
                nbinsx: 30
            }};
            Plotly.newPlot('cv', [trace1, trace2], {{
                ...darkTemplate.layout,
                title: 'Coefficient of Variation (CV) Distribution',
                xaxis: {{ title: 'CV' }},
                yaxis: {{ title: 'Frequency' }},
                barmode: 'overlay'
            }});
        }}

        function drawHeatmap() {{
            if (!plotData.heatmap_data) {{
                document.getElementById('heatmap').innerHTML = '<p>Not enough data for heatmap.</p>';
                return;
            }}
            const hm = plotData.heatmap_data;
            const trace = {{
                z: hm.values,
                x: hm.columns,
                y: hm.genes,
                type: 'heatmap',
                colorscale: 'Viridis',
                colorbar: {{ title: 'Rel. Usage' }}
            }};
            Plotly.newPlot('heatmap', [trace], {{
                ...darkTemplate.layout,
                title: 'APA Usage Heatmap (Clustered)',
                height: Math.max(600, hm.genes.length * 10),
                xaxis: {{ title: 'APA Site' }},
                yaxis: {{ title: 'Gene', automargin: true }}
            }});
        }}

        function drawTSNE() {{
            if (!plotData.tsne) {{
                document.getElementById('tsne').innerHTML = '<p>t-SNE not available (sklearn missing or insufficient genes).</p>';
                return;
            }}
            const tsne = plotData.tsne;
            const trace = {{
                x: tsne.x,
                y: tsne.y,
                mode: 'markers',
                type: 'scatter',
                text: tsne.genes,
                marker: {{
                    size: 8,
                    color: tsne.clusters,
                    colorscale: 'Portland',
                    showscale: true,
                    colorbar: {{ title: 'Cluster' }}
                }}
            }};
            Plotly.newPlot('tsne', [trace], {{
                ...darkTemplate.layout,
                title: 't-SNE Projection of Genes (Colored by Cluster)',
                xaxis: {{ title: 't-SNE 1' }},
                yaxis: {{ title: 't-SNE 2' }}
            }});
        }}

        // ---------- 初始化 ----------
        window.onload = function() {{
            drawChromosomeCharts();
            populateChromSelect();
            drawUsageHistogram('all');
            document.getElementById('chrom-select').addEventListener('change', e => drawUsageHistogram(e.target.value));

            drawSiteCountDist();
            drawMajorUsageDist();
            drawCVDist();
            drawHeatmap();
            drawTSNE();
        }};
    </script>
</body>
</html>
"""
    report_path = os.path.join(outdir, f"{prefix}_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"HTML 报告已生成: {report_path}")


def create_sqlite_db(outdir, prefix):
    """将关键表格存入 SQLite"""
    db_path = os.path.join(outdir, f"{prefix}_data.db")
    conn = sqlite3.connect(db_path)
    for fname in [f"{prefix}_apa_gene_level_summary.txt", f"{prefix}_apa_detailed_analysis.txt",
                  f"{prefix}_summary.txt", f"{prefix}_apa_summary.txt"]:
        path = os.path.join(outdir, fname)
        if os.path.exists(path):
            df = pd.read_csv(path, sep='\t')
            table_name = fname.replace('.txt', '').replace(f"{prefix}_", '')
            df.to_sql(table_name, conn, if_exists='replace', index=False)
    conn.close()
    print(f"SQLite 数据库已生成: {db_path}")


# ---------- 主函数 ----------
def main():
    parser = argparse.ArgumentParser(description='APA 分析完整流水线（全前端图表版）')
    parser.add_argument('--all_results', required=True, help='all_apa_results.txt 路径')
    parser.add_argument('--predictions', required=True, help='apa_predictions.txt 路径')
    parser.add_argument('--utr', required=True, help='refined_3utr.bed 路径')
    parser.add_argument('--genome', required=True, help='参考基因组 FASTA 文件路径')
    parser.add_argument('--outdir', default='./apa_output', help='输出目录')
    parser.add_argument('--prefix', default='apa', help='输出文件前缀')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. apa_summary
    apa_summary(args.all_results, args.outdir, args.prefix)

    # 2. deapa_module
    deapa_module(args.predictions, args.utr, args.genome, args.outdir, args.prefix)

    # 3. further_module (生成 JSON 数据)
    further_module(args.all_results, args.outdir, args.prefix)

    # 4. 生成 HTML 报告
    plot_json = os.path.join(args.outdir, f"{args.prefix}_plot_data.json")
    generate_html_report(args.outdir, args.prefix, args.predictions, args.utr, plot_json)

    # 5. 生成 SQLite 数据库
    create_sqlite_db(args.outdir, args.prefix)

    print("\n所有分析完成！请打开以下文件查看交互式报告：")
    print(f"   {os.path.join(args.outdir, args.prefix + '_report.html')}")


if __name__ == '__main__':
    main()