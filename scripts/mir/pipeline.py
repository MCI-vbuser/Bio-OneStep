#!/usr/bin/env python3
"""
pipeline.py – miRNA target energy analysis and extreme value fitting (memory-safe)
"""
import argparse, re, json, sqlite3
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gumbel_l, weibull_min, genextreme
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

MAX_SAMPLE_FOR_FIT = 10000  # 拟合时最大采样点，避免计算过慢

def extract_energies(input_path):
    """流式解析 miranda 输出，返回能量列表"""
    energies = []
    header_found = False
    idx = None
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # 寻找表头（可能出现在任意位置）
            if not header_found and ('Tot Energy' in line or 'Tot_Energy' in line):
                cols = re.split(r'[,\t]+', line)
                for candidate in ['Tot Energy', 'Tot_Energy']:
                    try:
                        idx = cols.index(candidate)
                        header_found = True
                        break
                    except ValueError:
                        continue
                if not header_found:
                    # 如果还是没找到，用默认索引 3
                    idx = 3
                    header_found = True
                continue
            # 提取 ">>" 行
            if line.startswith('>>'):
                parts = line[2:].split()
                if idx is not None and idx < len(parts):
                    try:
                        energy = float(parts[idx])
                        energies.append(energy)
                    except ValueError:
                        pass
    return energies

def neg_log_likelihood(params, dist, data, trunc_upper):
    # 保持不变
    if dist == 'gumbel_l':
        loc, scale = params
        if scale <= 0: return np.inf
        cdf_trunc = gumbel_l.cdf(trunc_upper, loc=loc, scale=scale)
        if cdf_trunc <= 0: return np.inf
        return -np.sum(gumbel_l.logpdf(data, loc=loc, scale=scale)) + len(data)*np.log(cdf_trunc)
    elif dist == 'weibull_min':
        c, loc, scale = params
        if scale <= 0 or c <= 0: return np.inf
        cdf_trunc = weibull_min.cdf(trunc_upper, c, loc=loc, scale=scale)
        if cdf_trunc <= 0: return np.inf
        return -np.sum(weibull_min.logpdf(data, c, loc=loc, scale=scale)) + len(data)*np.log(cdf_trunc)
    elif dist == 'genextreme':
        c, loc, scale = params
        if scale <= 0: return np.inf
        cdf_trunc = genextreme.cdf(trunc_upper, c, loc=loc, scale=scale)
        if cdf_trunc <= 0: return np.inf
        return -np.sum(genextreme.logpdf(data, c, loc=loc, scale=scale)) + len(data)*np.log(cdf_trunc)

def fit_distribution(dist_name, init, bounds, data, trunc_upper):
    res = minimize(neg_log_likelihood, x0=init, args=(dist_name, data, trunc_upper),
                   method='L-BFGS-B', bounds=bounds)
    if not res.success:
        return None, None, None
    params = res.x
    n = len(data)
    sorted_data = np.sort(data)
    p = (np.arange(1, n+1) - 0.5) / n
    # 计算理论分位数
    if dist_name == 'gumbel_l':
        loc, scale = params
        cdf_thr = gumbel_l.cdf(trunc_upper, loc=loc, scale=scale)
        theo = gumbel_l.ppf(p * cdf_thr, loc=loc, scale=scale)
    elif dist_name == 'weibull_min':
        c, loc, scale = params
        cdf_thr = weibull_min.cdf(trunc_upper, c, loc=loc, scale=scale)
        theo = weibull_min.ppf(p * cdf_thr, c, loc=loc, scale=scale)
    elif dist_name == 'genextreme':
        c, loc, scale = params
        cdf_thr = genextreme.cdf(trunc_upper, c, loc=loc, scale=scale)
        theo = genextreme.ppf(p * cdf_thr, c, loc=loc, scale=scale)
    else:
        return None, None, None
    r2 = np.corrcoef(sorted_data, theo)[0, 1] ** 2
    return params, r2, theo

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()

    print("流式提取自由能数据...")
    energies = extract_energies(args.input)
    print(f"共提取 {len(energies)} 条记录")
    if not energies:
        print("未提取到能量数据，退出。")
        return

    energy = np.array(energies)
    threshold = -5.0
    lower_bound = -100.0
    filtered = energy[(energy >= lower_bound) & (energy <= threshold)]
    print(f"过滤后数据量: {len(filtered)}")

    # 存入数据库（完整数据）
    db_path = f"{args.output_dir}/mirna_results.db"
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS energy_data")
    conn.execute("CREATE TABLE energy_data (Tot_Energy REAL)")
    for i in range(0, len(filtered), 10000):   # 分批写入，避免大事务
        batch = filtered[i:i+10000]
        conn.executemany("INSERT INTO energy_data VALUES (?)", [(e,) for e in batch])
    conn.commit()

    # 拟合时若数据量过大，随机采样
    if len(filtered) > MAX_SAMPLE_FOR_FIT:
        rng = np.random.RandomState(42)
        indices = rng.choice(len(filtered), MAX_SAMPLE_FOR_FIT, replace=False)
        fit_data = filtered[indices]
        print(f"从 {len(filtered)} 条中采样 {MAX_SAMPLE_FOR_FIT} 条用于拟合。")
    else:
        fit_data = filtered

    # 初始参数估计
    loc0_g, scale0_g = gumbel_l.fit(fit_data)
    try:
        c0_w, loc0_w, scale0_w = weibull_min.fit(fit_data)
    except:
        loc0_w = fit_data.min() - 0.1
        p63 = np.percentile(fit_data, 63.2)
        scale0_w = p63 - loc0_w
        c0_w = 1.5
    c0_gev, loc0_gev, scale0_gev = genextreme.fit(fit_data)

    configs = [
        ('Gumbel (left)', 'gumbel_l', [loc0_g, scale0_g], [(None,None),(1e-6,None)]),
        ('Weibull (min)', 'weibull_min', [c0_w,loc0_w,scale0_w], [(1e-6,None),(None,None),(1e-6,None)]),
        ('GEV', 'genextreme', [c0_gev,loc0_gev,scale0_gev], [(None,None),(None,None),(1e-6,None)])
    ]
    results = {}
    for name, dist, init, bounds in configs:
        print(f"拟合 {name} ...")
        params, r2, theo = fit_distribution(dist, init, bounds, fit_data, threshold)
        if params is not None:
            results[name] = {'params': params.tolist(), 'r2': r2, 'theoretical': theo}
            print(f"  R² = {r2:.4f}, 参数 = {params}")

    if not results:
        print("所有分布拟合失败。")
        conn.close()
        return

    best_name = max(results, key=lambda k: results[k]['r2'])
    best = results[best_name]
    summary = {
        'best_distribution': best_name,
        'R_squared': best['r2'],
        'parameters': best['params'],
        'data_count': len(filtered),
        'energy_range': [float(filtered.min()), float(filtered.max())]
    }
    json_path = f"{args.output_dir}/fit_summary.json"
    with open(json_path, 'w') as jf:
        json.dump(summary, jf, indent=2)

    # 写入拟合摘要表
    conn.execute("DROP TABLE IF EXISTS fit_summary")
    conn.execute("CREATE TABLE fit_summary (distribution TEXT, params TEXT, r2 REAL)")
    for name, res in results.items():
        conn.execute("INSERT INTO fit_summary VALUES (?,?,?)", (name, json.dumps(res['params']), res['r2']))
    conn.commit()
    conn.close()

    # 绘图（也基于采样数据，避免太密）
    fig, axes = plt.subplots(1, 2, figsize=(14,5))
    ax = axes[0]
    bin_width = 2.5
    n_bins = int((threshold - lower_bound) / bin_width)   # 结果 = 38
    ax.hist(fit_data, bins=n_bins, density=True, alpha=0.6, color='gray', edgecolor='black', label='Observed')
    x_grid = np.linspace(lower_bound, threshold, 500)
    colors = {'Gumbel (left)':'red','Weibull (min)':'blue','GEV':'green'}
    for name, res in results.items():
        params = res['params']
        if name == 'Gumbel (left)':
            loc,scale = params
            pdf = gumbel_l.pdf(x_grid, loc=loc, scale=scale)
            cdf_tr = gumbel_l.cdf(threshold, loc=loc, scale=scale)
        elif name == 'Weibull (min)':
            c,loc,scale = params
            pdf = weibull_min.pdf(x_grid, c, loc=loc, scale=scale)
            cdf_tr = weibull_min.cdf(threshold, c, loc=loc, scale=scale)
        else:
            c,loc,scale = params
            pdf = genextreme.pdf(x_grid, c, loc=loc, scale=scale)
            cdf_tr = genextreme.cdf(threshold, c, loc=loc, scale=scale)
        pdf_adj = np.where(x_grid <= threshold, pdf / cdf_tr, np.nan)
        ax.plot(x_grid, pdf_adj, color=colors[name], lw=2, label=f"{name} (R²={res['r2']:.3f})")
    ax.axvline(threshold, color='black', linestyle='--', label=f'Threshold {threshold}')
    ax.set_xlim(lower_bound, threshold)
    ax.legend(); ax.set_xlabel('Free Energy'); ax.set_ylabel('Density')
    ax.set_title('Truncated Extreme Value Fitting')

    ax2 = axes[1]
    sorted_data = np.sort(fit_data)
    theo_best = best['theoretical']
    ax2.scatter(theo_best, sorted_data, alpha=0.4, s=5, color='steelblue')
    mi, ma = min(theo_best.min(), sorted_data.min()), max(theo_best.max(), sorted_data.max())
    ax2.plot([mi, ma], [mi, ma], 'k--')
    ax2.set_xlabel('Theoretical Quantiles'); ax2.set_ylabel('Observed Quantiles')
    ax2.set_title(f'Q-Q Plot ({best_name}, R²={best["r2"]:.3f})')
    ax2.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    img_path = f"{args.output_dir}/energy_fit.png"
    plt.savefig(img_path, dpi=150)
    plt.close()
    print(f"拟合图已保存: {img_path}")
    generate_html_report(args.output_dir, db_path, summary)

def generate_html_report(output_dir, db_path, fit_summary):
    """
    生成交互式 HTML 报告：
    - 读取数据库中的 Tot_Energy 数据
    - 计算直方图（bins 和 counts）
    - 嵌入拟合结果（分布名称、R²、参数）
    - 使用 Chart.js 绘制直方图与拟合密度曲线的叠加图
    """
    import json
    conn = sqlite3.connect(db_path)
    cursor = conn.execute("SELECT Tot_Energy FROM energy_data")
    energies = np.array([row[0] for row in cursor], dtype=np.float64)
    conn.close()

    # 直方图设置
    bin_edges = np.arange(-100, -4, 2.5)  # [-80,-75), ..., [-10,-5]
    counts, _ = np.histogram(energies, bins=bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # 拟合密度曲线数据（用于 Chart.js）
    x_min, x_max = -100, 10
    x_grid = np.linspace(x_min, x_max, 300)
    best_dist = fit_summary['best_distribution']
    params = fit_summary['parameters']

    # 根据最佳分布计算 pdf
    if best_dist == 'Gumbel (left)':
        loc, scale = params
        pdf = gumbel_l.pdf(x_grid, loc=loc, scale=scale)
    elif best_dist == 'Weibull (min)':
        c, loc, scale = params
        pdf = weibull_min.pdf(x_grid, c, loc=loc, scale=scale)
    elif best_dist == 'GEV':
        c, loc, scale = params
        pdf = genextreme.pdf(x_grid, c, loc=loc, scale=scale)
    else:
        pdf = np.zeros_like(x_grid)

    # 数据点转成列表，准备嵌入 JS
    hist_data = {
        'centers': bin_centers.tolist(),
        'counts': counts.tolist()
    }
    curve_data = {
        'x': x_grid.tolist(),
        'y': pdf.tolist()
    }

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Miranda Energy Distribution Fitting</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <style>
        body {{ background: #0d1117; color: #c9d1d9; font-family: sans-serif; padding: 20px; }}
        h1, h2 {{ color: #f0f6fc; }}
        .card {{ background: #161b22; border: 1px solid #30363d; border-radius: 6px; padding: 16px; margin: 12px; display: inline-block; }}
        .card h3 {{ color: #8b949e; margin: 0 0 8px; }}
        .card .value {{ font-size: 2em; font-weight: bold; }}
        .chart-container {{ height: 400px; margin: 20px 0; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #30363d; padding: 8px; text-align: left; }}
        th {{ background: #21262d; }}
    </style>
</head>
<body>
    <h1>miRNA Target Energy Extreme Value Fitting</h1>
    <div class="card"><h3>Best Distribution</h3><div class="value">{fit_summary['best_distribution']}</div></div>
    <div class="card"><h3>R²</h3><div class="value">{fit_summary['R_squared']:.4f}</div></div>
    <div class="card"><h3>Parameters</h3><div class="value">{', '.join(f'{p:.4f}' for p in fit_summary['parameters'])}</div></div>
    <div class="card"><h3>Data Points</h3><div class="value">{fit_summary['data_count']}</div></div>

    <h2>Energy Histogram & Fitted Density</h2>
    <div class="chart-container"><canvas id="mainChart"></canvas></div>

    <h2>Fit Parameters Details</h2>
    <table>
        <tr><th>Distribution</th><th>Parameters</th><th>R²</th></tr>
        {''.join(f"<tr><td>{name}</td><td>{', '.join(f'{v:.4f}' for v in res['params'])}</td><td>{res['r2']:.4f}</td></tr>"
                 for name, res in [('Best', {'params': params, 'r2': fit_summary['R_squared']})])}
    </table>

    <script>
        const histData = {json.dumps(hist_data)};
        const curveData = {json.dumps(curve_data)};
        const threshold = -5.0;  // -5.0

        const ctx = document.getElementById('mainChart').getContext('2d');
        new Chart(ctx, {{
            type: 'bar',
            data: {{
                labels: histData.centers.map(c => `${{c.toFixed(0)}}`),
                datasets: [
                    {{
                        label: 'Frequency',
                        data: histData.counts,
                        backgroundColor: 'rgba(128,128,128,0.6)',
                        borderColor: '#808080',
                        borderWidth: 1,
                        yAxisID: 'y'
                    }},
                    {{
                        label: `Fitted ${{'{fit_summary["best_distribution"]}'}}`,
                        data: curveData.y,
                        backgroundColor: 'transparent',
                        borderColor: '#ff6384',
                        borderWidth: 2,
                        type: 'line',
                        yAxisID: 'y1',
                        pointRadius: 0,
                        tension: 0.2,
                        spanGaps: true
                    }}
                ]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    x: {{ title: {{ display: true, text: 'Tot Energy (kcal/mol)', color: '#c9d1d9' }}, ticks: {{ color: '#c9d1d9' }} }},
                    y: {{ title: {{ display: true, text: 'Count', color: '#c9d1d9' }}, position: 'left', ticks: {{ color: '#c9d1d9' }}, grid: {{ color: '#30363d' }} }},
                    y1: {{ title: {{ display: true, text: 'Density', color: '#c9d1d9' }}, position: 'right', ticks: {{ color: '#c9d1d9' }}, grid: {{ display: false }} }}
                }},
                plugins: {{
                    legend: {{ labels: {{ color: '#c9d1d9' }} }},
                    tooltip: {{ mode: 'index', intersect: false }}
                }}
            }}
        }});
    </script>
</body>
</html>"""

    report_path = os.path.join(output_dir, "mirna_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"HTML 报告已保存：{report_path}")

if __name__ == '__main__':
    main()