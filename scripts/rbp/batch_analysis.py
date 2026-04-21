#!/usr/bin/env python
"""
批量分析脚本 - 处理多个样本和多个RBPs
"""
import os
import glob
import pandas as pd
from rbp_inference import RBPInference

def analyze_sample_files(sample_dir, output_dir):
    """
    批量分析样本文件
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 初始化推理器
    inferencer = RBPInference()
    
    # 查找样本文件
    fasta_files = glob.glob(os.path.join(sample_dir, "*.fa")) + \
                  glob.glob(os.path.join(sample_dir, "*.fasta"))
    bed_files = glob.glob(os.path.join(sample_dir, "*.bed"))
    
    all_results = []
    
    # 处理FASTA文件
    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file}")
        sample_name = os.path.basename(fasta_file).split('.')[0]
        
        results = inferencer.analyze_novel_transcripts(fasta_file)
        results['sample'] = sample_name
        
        # 保存单独结果
        output_file = os.path.join(output_dir, f"{sample_name}_rbp_predictions.csv")
        results.to_csv(output_file, index=False)
        
        all_results.append(results)
    
    # 合并所有结果
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        combined_results.to_csv(
            os.path.join(output_dir, "all_samples_rbp_predictions.csv"), 
            index=False
        )
        
        # 生成摘要报告
        generate_summary_report(combined_results, output_dir)
        
        return combined_results
    else:
        print("No valid FASTA files found in sample directory.")
        return pd.DataFrame()

def generate_summary_report(results_df, output_dir):
    """
    生成分析摘要报告
    """
    if results_df.empty:
        print("No results to generate summary report.")
        return
    
    # 按RBP统计
    rbp_stats = results_df.groupby('rbp').agg({
        'score': ['mean', 'max', 'min'],
        'binding': 'sum',
        'sequence': 'count'
    }).round(3)
    
    rbp_stats.columns = ['avg_score', 'max_score', 'min_score', 
                        'binding_count', 'total_sequences']
    rbp_stats['binding_percent'] = (rbp_stats['binding_count'] / 
                                   rbp_stats['total_sequences'] * 100).round(1)
    
    rbp_stats = rbp_stats.sort_values('avg_score', ascending=False)
    
    # 保存统计结果
    stats_file = os.path.join(output_dir, "rbp_summary_statistics.csv")
    rbp_stats.to_csv(stats_file)
    
    # 生成热力图数据
    pivot_table = results_df.pivot_table(
        index='sequence',
        columns='rbp',
        values='score',
        aggfunc='mean'
    ).fillna(0)
    
    heatmap_file = os.path.join(output_dir, "rbp_binding_heatmap.csv")
    pivot_table.to_csv(heatmap_file)
    
    print(f"Summary statistics saved to: {stats_file}")
    print(f"Binding heatmap data saved to: {heatmap_file}")
    
    return rbp_stats

def analyze_with_gtf(gtf_file, genome_fasta, output_dir, gene_list=None):
    """
    使用GTF文件分析特定基因的UTR区域
    """
    from pybedtools import BedTool
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取GTF文件
    gtf = BedTool(gtf_file)
    
    # 筛选3'UTR区域
    utr_regions = gtf.filter(lambda x: x[2] == 'three_prime_utr')
    
    if gene_list:
        # 筛选特定基因
        gene_bed = []
        for region in utr_regions:
            gene_name = region.attrs.get('gene_name', '')
            if gene_name in gene_list:
                gene_bed.append(region)
        utr_regions = BedTool(gene_bed)
    
    # 保存为临时BED文件
    temp_bed = os.path.join(output_dir, "temp_utr_regions.bed")
    utr_regions.saveas(temp_bed)
    
    # 进行RBP预测
    inferencer = RBPInference()
    results = inferencer.predict_bed_regions(temp_bed, genome_fasta)
    
    if results is not None:
        output_file = os.path.join(output_dir, "utr_rbp_predictions.csv")
        results.to_csv(output_file, index=False)
        
        # 清理临时文件
        os.remove(temp_bed)
        
        print(f"UTR analysis completed. Results saved to: {output_file}")
        
        return results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='批量RBP结合分析')
    parser.add_argument('--input-dir', type=str, required=True,
                       help='输入文件目录')
    parser.add_argument('--output-dir', type=str, default='./rbp_results',
                       help='输出目录')
    parser.add_argument('--gtf', type=str,
                       help='GTF文件路径（可选）')
    parser.add_argument('--genome', type=str,
                       help='基因组FASTA文件路径')
    parser.add_argument('--genes', nargs='+',
                       help='要分析的基因列表')
    
    args = parser.parse_args()
    
    # 根据输入类型选择分析方法
    if args.gtf and args.genome:
        print(f"Analyzing UTR regions from GTF: {args.gtf}")
        results = analyze_with_gtf(
            args.gtf, 
            args.genome, 
            args.output_dir,
            args.genes
        )
    else:
        print(f"Batch analyzing files in: {args.input_dir}")
        results = analyze_sample_files(args.input_dir, args.output_dir)
