#!/usr/bin/env python
"""
RBP分析结果可视化
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import matplotlib as mpl
import os

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

def plot_rbp_binding_frequency(results_df, top_n=20, save_path=None):
    """
    绘制RBP结合频率条形图
    """
    # 计算每个RBP的结合频率
    binding_stats = results_df.groupby('rbp')['binding'].agg([
        'sum', 'count'
    ]).reset_index()
    binding_stats['frequency'] = binding_stats['sum'] / binding_stats['count'] * 100
    
    # 按频率排序并选择top_n
    binding_stats = binding_stats.sort_values('frequency', ascending=False).head(top_n)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 创建颜色映射
    colors = cm.viridis(np.linspace(0.2, 0.8, len(binding_stats)))
    
    # 绘制条形图
    bars = ax.barh(range(len(binding_stats)), 
                   binding_stats['frequency'], 
                   color=colors)
    
    # 添加数值标签
    for i, (bar, row) in enumerate(zip(bars, binding_stats.itertuples())):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{row.frequency:.1f}% (n={row.sum}/{row.count})',
                va='center', fontsize=9)
    
    # 设置y轴标签
    ax.set_yticks(range(len(binding_stats)))
    ax.set_yticklabels(binding_stats['rbp'], fontsize=10)
    
    # 设置标题和标签
    ax.set_xlabel('Binding Frequency (%)', fontsize=12)
    ax.set_title(f'Top {top_n} RBP Binding Frequencies', fontsize=14, fontweight='bold')
    
    # 调整布局
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Binding frequency plot saved to: {save_path}")
    
    plt.show()
    return fig, ax

def plot_score_distribution(results_df, rbp_list=None, save_path=None):
    """
    绘制多个RBP的分数分布图
    """
    if rbp_list is None:
        # 选择分数最高的前10个RBP
        avg_scores = results_df.groupby('rbp')['score'].mean().sort_values(ascending=False)
        rbp_list = avg_scores.head(10).index.tolist()
    
    # 筛选数据
    filtered_data = results_df[results_df['rbp'].isin(rbp_list)]
    
    # 创建图表
    fig, axes = plt.subplots(2, 5, figsize=(20, 10))
    axes = axes.flatten()
    
    for idx, rbp in enumerate(rbp_list):
        if idx >= len(axes):
            break
            
        rbp_data = filtered_data[filtered_data['rbp'] == rbp]
        
        # 分离正负样本
        positive_scores = rbp_data[rbp_data['binding']]['score']
        negative_scores = rbp_data[~rbp_data['binding']]['score']
        
        # 绘制分布
        ax = axes[idx]
        if len(positive_scores) > 0:
            sns.kdeplot(positive_scores, ax=ax, label='Positive', color='red', fill=True, alpha=0.5)
        if len(negative_scores) > 0:
            sns.kdeplot(negative_scores, ax=ax, label='Negative', color='blue', fill=True, alpha=0.5)
        
        ax.set_title(f'{rbp}\n(n={len(rbp_data)})', fontsize=10, fontweight='bold')
        ax.set_xlabel('Score')
        ax.set_xlim(0, 1)
        ax.legend(fontsize=8)
    
    # 隐藏多余的子图
    for idx in range(len(rbp_list), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Score Distribution for Top RBPs', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Score distribution plot saved to: {save_path}")
    
    plt.show()
    return fig, axes

def plot_heatmap(results_df, save_path=None):
    """
    绘制RBP结合热力图
    """
    # 创建数据透视表
    pivot_data = results_df.pivot_table(
        index='sequence',
        columns='rbp',
        values='score',
        aggfunc='mean'
    ).fillna(0)
    
    # 只选择有数据的RBP和序列
    pivot_data = pivot_data.loc[pivot_data.sum(axis=1) > 0, 
                               pivot_data.sum(axis=0) > 0]
    
    # 如果数据太多，选择前50行和前30列
    if pivot_data.shape[0] > 50:
        pivot_data = pivot_data.head(50)
    if pivot_data.shape[1] > 30:
        pivot_data = pivot_data[pivot_data.columns[:30]]
    
    # 创建热力图
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # 使用seaborn绘制热力图
    sns.heatmap(pivot_data, 
                cmap='YlOrRd',
                ax=ax,
                cbar_kws={'label': 'Binding Score'},
                square=True)
    
    ax.set_title('RBP Binding Score Heatmap', fontsize=16, fontweight='bold')
    ax.set_xlabel('RNA Binding Proteins', fontsize=12)
    ax.set_ylabel('Sequences/Regions', fontsize=12)
    
    # 调整x轴标签角度
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved to: {save_path}")
    
    plt.show()
    return fig, ax

def create_comprehensive_report(results_dir, output_dir='./reports'):
    """
    创建综合分析报告
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 查找所有结果文件
    result_files = glob.glob(os.path.join(results_dir, "*.csv"))
    
    all_results = []
    for file in result_files:
        try:
            df = pd.read_csv(file)
            df['source_file'] = os.path.basename(file)
            all_results.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if not all_results:
        print("No valid result files found.")
        return
    
    # 合并所有结果
    combined_results = pd.concat(all_results, ignore_index=True)
    
    # 创建报告目录
    report_dir = os.path.join(output_dir, 'figures')
    os.makedirs(report_dir, exist_ok=True)
    
    # 生成可视化
    print("Generating visualizations...")
    
    # 1. RBP结合频率图
    freq_fig, _ = plot_rbp_binding_frequency(
        combined_results,
        top_n=25,
        save_path=os.path.join(report_dir, 'rbp_binding_frequency.png')
    )
    
    # 2. 分数分布图
    score_fig, _ = plot_score_distribution(
        combined_results,
        save_path=os.path.join(report_dir, 'score_distribution.png')
    )
    
    # 3. 热力图
    heatmap_fig, _ = plot_heatmap(
        combined_results,
        save_path=os.path.join(report_dir, 'binding_heatmap.png')
    )
    
    # 生成文本报告
    generate_text_report(combined_results, output_dir)
    
    print(f"Comprehensive report generated in: {output_dir}")
    return combined_results

def generate_text_report(results_df, output_dir):
    """
    生成文本格式的分析报告
    """
    report_file = os.path.join(output_dir, 'analysis_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("RBP BINDING ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        # 基本信息
        f.write("1. BASIC INFORMATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total predictions: {len(results_df):,}\n")
        f.write(f"Unique RBPs: {results_df['rbp'].nunique()}\n")
        f.write(f"Unique sequences: {results_df['sequence'].nunique()}\n")
        f.write(f"Model types: {results_df['model_type'].unique()}\n")
        f.write(f"Samples: {results_df.get('sample', 'N/A')}\n\n")
        
        # 结合统计
        f.write("2. BINDING STATISTICS\n")
        f.write("-" * 40 + "\n")
        total_binding = results_df['binding'].sum()
        binding_percent = total_binding / len(results_df) * 100
        f.write(f"Total binding events: {total_binding:,} ({binding_percent:.1f}%)\n")
        f.write(f"Non-binding events: {len(results_df) - total_binding:,} ({100 - binding_percent:.1f}%)\n\n")
        
        # Top RBPs
        f.write("3. TOP 15 RBPs BY BINDING FREQUENCY\n")
        f.write("-" * 40 + "\n")
        
        binding_stats = results_df.groupby('rbp').agg({
            'binding': ['sum', 'count'],
            'score': 'mean'
        }).round(3)
        
        binding_stats.columns = ['binding_count', 'total_count', 'avg_score']
        binding_stats['frequency'] = (binding_stats['binding_count'] / 
                                     binding_stats['total_count'] * 100).round(1)
        
        top_rbps = binding_stats.sort_values('frequency', ascending=False).head(15)
        
        f.write(f"{'RBP':<15} {'Binding':>10} {'Total':>10} {'Freq%':>8} {'AvgScore':>10}\n")
        f.write("-" * 60 + "\n")
        
        for idx, row in top_rbps.iterrows():
            f.write(f"{idx:<15} {row['binding_count']:>10,} {row['total_count']:>10,} "
                   f"{row['frequency']:>7.1f}% {row['avg_score']:>9.3f}\n")
        
        f.write("\n")
        
        # 按模型类型统计
        f.write("4. STATISTICS BY MODEL TYPE\n")
        f.write("-" * 40 + "\n")
        
        for model_type in ['high', 'med', 'low']:
            if model_type in results_df['model_type'].unique():
                model_data = results_df[results_df['model_type'] == model_type]
                model_binding = model_data['binding'].sum()
                model_percent = model_binding / len(model_data) * 100
                
                f.write(f"{model_type.upper()} model:\n")
                f.write(f"  Predictions: {len(model_data):,}\n")
                f.write(f"  Binding events: {model_binding:,} ({model_percent:.1f}%)\n")
                f.write(f"  RBPs in model: {model_data['rbp'].nunique()}\n")
        
        f.write("\n")
        
        # 高置信度预测
        f.write("5. HIGH-CONFIDENCE PREDICTIONS (Score > 0.9)\n")
        f.write("-" * 40 + "\n")
        
        high_conf = results_df[results_df['score'] > 0.9]
        if len(high_conf) > 0:
            f.write(f"High-confidence predictions: {len(high_conf):,}\n")
            
            # 高置信度的Top RBPs
            high_conf_stats = high_conf.groupby('rbp').agg({
                'sequence': 'count',
                'score': 'mean'
            }).sort_values('sequence', ascending=False).head(10)
            
            for rbp, row in high_conf_stats.iterrows():
                f.write(f"  {rbp:<15}: {row['sequence']:>4,} predictions, avg score: {row['score']:.3f}\n")
        else:
            f.write("No high-confidence predictions found.\n")
    
    print(f"Text report saved to: {report_file}")

if __name__ == "__main__":
    import argparse
    import glob
    
    parser = argparse.ArgumentParser(description='RBP分析结果可视化')
    parser.add_argument('--input', type=str, required=True,
                       help='输入CSV文件或包含CSV文件的目录')
    parser.add_argument('--output', type=str, default='./reports',
                       help='输出目录')
    parser.add_argument('--create-report', action='store_true',
                       help='创建综合分析报告')
    
    args = parser.parse_args()
    
    if os.path.isdir(args.input):
        create_comprehensive_report(args.input, args.output)
    elif os.path.isfile(args.input):
        # 加载单个文件
        results_df = pd.read_csv(args.input)
        
        # 创建输出目录
        output_dir = args.output
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成可视化
        plot_rbp_binding_frequency(
            results_df,
            save_path=os.path.join(output_dir, 'binding_frequency.png')
        )
        
        plot_score_distribution(
            results_df,
            save_path=os.path.join(output_dir, 'score_distribution.png')
        )
        
        if len(results_df) < 100000:  # 避免太大的热力图
            plot_heatmap(
                results_df,
                save_path=os.path.join(output_dir, 'heatmap.png')
            )
        
        generate_text_report(results_df, output_dir)
    else:
        print(f"Error: Input path {args.input} not found.")
