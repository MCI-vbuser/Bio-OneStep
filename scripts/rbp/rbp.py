#!/usr/bin/env python3
"""
RBP 分析报告生成器（数据库 + 交互式图表）
用法: python rbp.py --input <csv_file> [--output <html_file>] [--db <db_file>]
"""

import argparse
import os
import sqlite3
import json
from datetime import datetime
import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

pio.templates.default = "plotly_dark"

# ---------- 数据库操作 ----------
def init_db(db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS rbp_predictions (
            id INTEGER PRIMARY KEY,
            sequence TEXT,
            rbp TEXT,
            model_type TEXT,
            score REAL,
            binding INTEGER
        )
    ''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS rbp_summary (
            rbp TEXT PRIMARY KEY,
            binding_count INTEGER,
            total_count INTEGER,
            avg_score REAL,
            frequency REAL
        )
    ''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS model_stats (
            model_type TEXT PRIMARY KEY,
            total_predictions INTEGER,
            binding_events INTEGER
        )
    ''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS co_binding (
            rbp1 TEXT,
            rbp2 TEXT,
            co_occurrence INTEGER,
            PRIMARY KEY (rbp1, rbp2)
        )
    ''')
    conn.commit()
    return conn

def save_to_db(conn, df):
    # 清空旧数据
    conn.execute("DELETE FROM rbp_predictions")
    conn.execute("DELETE FROM rbp_summary")
    conn.execute("DELETE FROM model_stats")
    conn.execute("DELETE FROM co_binding")

    # 确保 binding 列为整数 (0/1)
    df_to_save = df.copy()
    df_to_save['binding'] = df_to_save['binding'].astype(int)
    df_to_save[['sequence', 'rbp', 'model_type', 'score', 'binding']].to_sql(
        'rbp_predictions', conn, if_exists='append', index=False
    )

    # 每个 RBP 的汇总统计
    rbp_stats = df.groupby('rbp').agg(
        binding_count=('binding', 'sum'),
        total_count=('binding', 'count'),
        avg_score=('score', 'mean')
    ).reset_index()
    rbp_stats['frequency'] = rbp_stats['binding_count'] / rbp_stats['total_count'] * 100
    rbp_stats.to_sql('rbp_summary', conn, if_exists='append', index=False)

    # 按模型类型统计
    model_stats = df.groupby('model_type').agg(
        total_predictions=('binding', 'count'),
        binding_events=('binding', 'sum')
    ).reset_index()
    model_stats.to_sql('model_stats', conn, if_exists='append', index=False)

    # 共结合分析（score > 0.7 且 binding=True）
    high_conf = df[(df['score'] > 0.7) & (df['binding'])]
    seq_rbps = defaultdict(set)
    for _, row in high_conf.iterrows():
        seq_rbps[row['sequence']].add(row['rbp'])

    co_occurrence = defaultdict(int)
    for rbps in seq_rbps.values():
        rbp_list = list(rbps)
        for i in range(len(rbp_list)):
            for j in range(i+1, len(rbp_list)):
                pair = tuple(sorted([rbp_list[i], rbp_list[j]]))
                co_occurrence[pair] += 1

    co_data = [(p[0], p[1], cnt) for p, cnt in co_occurrence.items()]
    if co_data:
        co_df = pd.DataFrame(co_data, columns=['rbp1', 'rbp2', 'co_occurrence'])
        co_df.to_sql('co_binding', conn, if_exists='append', index=False)

    conn.commit()

# ---------- 图表函数 ----------
def plot_binding_frequency(rbp_stats):
    top = rbp_stats.sort_values('frequency', ascending=False).head(20)
    fig = go.Figure(data=[
        go.Bar(
            x=top['frequency'],
            y=top['rbp'],
            orientation='h',
            text=top['frequency'].apply(lambda x: f'{x:.1f}%'),
            textposition='outside',
            marker=dict(color=top['frequency'], colorscale='Viridis'),
            hovertemplate='<b>%{y}</b><br>Frequency: %{x:.1f}%<br>Binding: %{customdata[0]}/%{customdata[1]}<extra></extra>',
            customdata=np.stack((top['binding_count'], top['total_count']), axis=-1)
        )
    ])
    fig.update_layout(
        title='RBP Binding Frequency (Top 20)',
        xaxis_title='Frequency (%)',
        yaxis_title='',
        height=600,
        margin=dict(l=150),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def plot_score_distribution(df, top_n=10):
    try:
        from scipy.stats import gaussian_kde
        use_kde = True
    except ImportError:
        use_kde = False

    top_rbps = df.groupby('rbp')['score'].mean().sort_values(ascending=False).head(top_n).index.tolist()
    df_sub = df[df['rbp'].isin(top_rbps)]

    rows = (top_n + 4) // 5
    cols = min(5, top_n)
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=top_rbps[:top_n],
                        shared_xaxes=True, shared_yaxes=True,
                        vertical_spacing=0.1, horizontal_spacing=0.1)

    for i, rbp in enumerate(top_rbps[:top_n]):
        row = i // cols + 1
        col = i % cols + 1
        data_rbp = df_sub[df_sub['rbp'] == rbp]

        for bind, color, name, fillcolor in [
            (True, 'red', 'positive', 'rgba(255,0,0,0.2)'),
            (False, 'blue', 'negative', 'rgba(0,0,255,0.2)')
        ]:
            scores = data_rbp[data_rbp['binding'] == bind]['score'].values
            if len(scores) < 2:
                continue
            if use_kde:
                x_grid = np.linspace(0, 1, 200)
                kde = gaussian_kde(scores)
                y_vals = kde.evaluate(x_grid)
                fig.add_trace(
                    go.Scatter(
                        x=x_grid, y=y_vals, mode='lines',
                        line=dict(color=color, width=2),
                        fill='tozeroy', fillcolor=fillcolor,
                        name=f'{rbp} {name}', showlegend=False
                    ),
                    row=row, col=col
                )
            else:
                fig.add_trace(
                    go.Histogram(
                        x=scores, histnorm='probability density',
                        marker_color=color, opacity=0.5,
                        name=f'{rbp} {name}', showlegend=False,
                        nbinsx=20
                    ),
                    row=row, col=col
                )
        fig.update_xaxes(range=[0, 1], title_text='Score' if row == rows else '', row=row, col=col)
        fig.update_yaxes(title_text='Density' if col == 1 else '', row=row, col=col)

    fig.update_layout(
        title=f'Score Distribution for Top {top_n} RBPs',
        height=300 * rows,
        showlegend=False,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def plot_heatmap(df):
    pivot = df.pivot_table(index='sequence', columns='rbp', values='score', aggfunc='mean').fillna(0)
    if pivot.shape[0] > 50:
        pivot = pivot.iloc[:50, :]
    if pivot.shape[1] > 30:
        pivot = pivot.iloc[:, :30]
    if pivot.empty:
        return None
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale='YlOrRd',
        colorbar=dict(title='Score')
    ))
    fig.update_layout(
        title='RBP Binding Score Heatmap',
        xaxis_title='RBP',
        yaxis_title='Sequence',
        height=600,
        xaxis=dict(tickangle=45),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def plot_network(df, min_co_occurrence=3):
    high_conf = df[(df['score'] > 0.7) & (df['binding'])]
    seq_rbps = defaultdict(set)
    for _, row in high_conf.iterrows():
        seq_rbps[row['sequence']].add(row['rbp'])

    co_occurrence = defaultdict(int)
    for rbps in seq_rbps.values():
        rbp_list = list(rbps)
        for i in range(len(rbp_list)):
            for j in range(i+1, len(rbp_list)):
                pair = tuple(sorted([rbp_list[i], rbp_list[j]]))
                co_occurrence[pair] += 1

    edges = [(p[0], p[1], w) for p, w in co_occurrence.items() if w >= min_co_occurrence]
    if not edges:
        return None

    G = nx.Graph()
    for rbp in df['rbp'].unique():
        G.add_node(rbp)
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)

    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    edge_trace = []
    for u, v, w in edges:
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_trace.append(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode='lines',
            line=dict(width=w/2, color='rgba(200,200,200,0.5)'),
            hoverinfo='none', showlegend=False
        ))

    binding_freq = df.groupby('rbp')['binding'].mean().to_dict()
    node_x, node_y, node_text, node_size, node_color = [], [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        freq = binding_freq.get(node, 0)
        node_size.append(20 + freq * 50)
        node_color.append(freq)
        node_text.append(f"{node}<br>Binding freq: {freq:.2f}")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[n for n in G.nodes()],
        textposition="top center",
        hovertext=node_text,
        marker=dict(size=node_size, color=node_color, colorscale='YlOrRd', showscale=True,
                    colorbar=dict(title='Binding Freq')),
        showlegend=False
    )
    fig = go.Figure(data=edge_trace + [node_trace])
    fig.update_layout(
        title='RBP Co-binding Network',
        showlegend=False, hovermode='closest',
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        height=700,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def generate_html_report(csv_path, html_path, db_path):
    print(f"读取数据: {csv_path}")
    df = pd.read_csv(csv_path)
    if 'binding' in df.columns:
        # 将 binding 列转换为布尔值（用于后续计算）
        df['binding'] = df['binding'].astype(bool)
    else:
        raise ValueError("CSV 缺少 'binding' 列")

    print(f"写入数据库: {db_path}")
    conn = init_db(db_path)
    save_to_db(conn, df)
    conn.close()

    total_predictions = len(df)
    unique_rbps = df['rbp'].nunique()
    unique_sequences = df['sequence'].nunique()
    total_binding = df['binding'].sum()
    binding_percent = total_binding / total_predictions * 100

    rbp_stats = df.groupby('rbp').agg(
        binding_count=('binding', 'sum'),
        total_count=('binding', 'count'),
        avg_score=('score', 'mean')
    ).reset_index()
    rbp_stats['frequency'] = rbp_stats['binding_count'] / rbp_stats['total_count'] * 100
    top_rbps = rbp_stats.sort_values('frequency', ascending=False).head(15)

    print("生成图表...")
    fig_freq = plot_binding_frequency(rbp_stats)
    fig_dist = plot_score_distribution(df)
    fig_heat = plot_heatmap(df)
    fig_net = plot_network(df)

    # 转换为 HTML div 片段
    div_freq = fig_freq.to_html(full_html=False, include_plotlyjs='cdn')
    div_dist = fig_dist.to_html(full_html=False, include_plotlyjs=False)
    div_heat = fig_heat.to_html(full_html=False, include_plotlyjs=False) if fig_heat else "<p>No heatmap data available.</p>"
    div_net = fig_net.to_html(full_html=False, include_plotlyjs=False) if fig_net else "<p>No network data (min co‑occurrence threshold not met).</p>"

    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>RBP Binding Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            background-color: #0d1117;
            color: #c9d1d9;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
            line-height: 1.5;
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        h1, h2, h3 {{
            color: #f0f6fc;
            border-bottom: 1px solid #30363d;
            padding-bottom: 0.3em;
        }}
        table.stat-table {{
            border-collapse: collapse;
            width: 100%;
            background-color: #161b22;
            border: 1px solid #30363d;
            border-radius: 6px;
            margin-bottom: 20px;
        }}
        table.stat-table td, table.stat-table th {{
            padding: 8px 12px;
            border: 1px solid #30363d;
        }}
        table.stat-table th {{
            background-color: #21262d;
            font-weight: 600;
        }}
        .plot-container {{
            margin: 30px 0;
            background-color: #161b22;
            border-radius: 6px;
            padding: 15px;
            border: 1px solid #30363d;
        }}
        .footer {{
            text-align: center;
            color: #8b949e;
            font-size: 0.9em;
            margin-top: 40px;
            border-top: 1px solid #30363d;
            padding-top: 20px;
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>RBP Binding Analysis Report</h1>

    <div class="plot-container">
        <h2>Summary Statistics</h2>
        <table class="stat-table">
            <tr><th>Total predictions</th><td>{total_predictions:,}</td></tr>
            <tr><th>Unique RBPs</th><td>{unique_rbps}</td></tr>
            <tr><th>Unique sequences</th><td>{unique_sequences:,}</td></tr>
            <tr><th>Binding events</th><td>{total_binding:,} ({binding_percent:.1f}%)</td></tr>
        </table>

        <h3>Top 15 RBPs by Binding Frequency</h3>
        <table class="stat-table">
            <tr><th>RBP</th><th>Binding</th><th>Total</th><th>Freq%</th><th>Avg Score</th></tr>
            {''.join(f'<tr><td>{row["rbp"]}</td><td>{row["binding_count"]:,}</td><td>{row["total_count"]:,}</td><td>{row["frequency"]:.1f}%</td><td>{row["avg_score"]:.3f}</td></tr>' for _, row in top_rbps.iterrows())}
        </table>
    </div>

    <div class="plot-container"><h2>Binding Frequency</h2>{div_freq}</div>
    <div class="plot-container"><h2>Score Distribution (Top 10 RBPs)</h2>{div_dist}</div>
    <div class="plot-container"><h2>Binding Score Heatmap</h2>{div_heat}</div>
    <div class="plot-container"><h2>RBP Co-binding Network</h2>{div_net}</div>

</div>
</body>
</html>"""

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print(f"HTML 报告已保存: {html_path}")
    print(f"数据库已保存: {db_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RBP analysis report with database and interactive plots")
    parser.add_argument('--input', required=True, help='Input CSV file (RBP predictions)')
    parser.add_argument('--output', default='./rbp_report.html', help='Output HTML file path')
    parser.add_argument('--db', default='./rbp/rbp_data.db', help='SQLite database file path')
    args = parser.parse_args()

    generate_html_report(args.input, args.output, args.db)