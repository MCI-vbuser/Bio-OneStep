#!/usr/bin/env python
"""
过滤RBP分析结果，保留binding为True的记录
支持大文件处理和多种过滤选项
"""
import pandas as pd
import numpy as np
import argparse
import os
import sys
from pathlib import Path
import gc

def filter_binding_true(input_file, output_file=None, 
                        min_score=0.0, top_n=None,
                        chunk_size=100000, verbose=True):
    """
    过滤CSV文件，保留binding为True的记录
    
    Parameters:
    -----------
    input_file : str
        输入CSV文件路径
    output_file : str
        输出CSV文件路径（如为None则自动生成）
    min_score : float
        最小分数阈值
    top_n : int
        只保留分数最高的前N条记录（按每个RBP）
    chunk_size : int
        分块处理大小（用于大文件）
    verbose : bool
        是否显示详细信息
    
    Returns:
    --------
    str : 输出文件路径
    int : 过滤后的记录数
    """
    # 自动生成输出文件名
    if output_file is None:
        input_path = Path(input_file)
        output_file = str(input_path.parent / f"{input_path.stem}_binding_true.csv")
    
    # 检查输入文件
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在: {input_file}")
        return None, 0
    
    if verbose:
        print(f"处理文件: {input_file}")
        print(f"输出文件: {output_file}")
    
    # 先读取前几行检查格式
    try:
        sample_df = pd.read_csv(input_file, nrows=5)
    except Exception as e:
        print(f"错误: 无法读取CSV文件: {e}")
        return None, 0
    
    # 检查必要的列
    required_cols = ['binding']
    missing_cols = [col for col in required_cols if col not in sample_df.columns]
    if missing_cols:
        print(f"错误: 缺少必要的列: {missing_cols}")
        print(f"可用列: {list(sample_df.columns)}")
        return None, 0
    
    # 判断文件大小，决定是否分块处理
    file_size_mb = os.path.getsize(input_file) / (1024 * 1024)
    use_chunks = file_size_mb > 50  # 大于50MB的文件使用分块处理
    
    if verbose:
        print(f"文件大小: {file_size_mb:.1f} MB")
        print(f"使用分块处理: {'是' if use_chunks else '否'}")
    
    total_rows = 0
    filtered_rows = 0
    
    if use_chunks:
        # 分块处理大文件
        print(f"使用分块处理，块大小: {chunk_size:,} 行")
        
        # 先读取所有块并收集binding为True的行
        chunks = []
        for i, chunk in enumerate(pd.read_csv(input_file, chunksize=chunk_size)):
            total_rows += len(chunk)
            
            # 过滤binding为True且分数大于阈值
            mask = chunk['binding'] == True
            if min_score > 0:
                mask = mask & (chunk['score'] >= min_score)
            
            filtered_chunk = chunk[mask].copy()
            
            if not filtered_chunk.empty:
                chunks.append(filtered_chunk)
                filtered_rows += len(filtered_chunk)
            
            if verbose and (i+1) % 10 == 0:
                print(f"  已处理 {total_rows:,} 行，保留 {filtered_rows:,} 行")
            
            # 释放内存
            del chunk
            gc.collect()
        
        # 合并所有块
        if chunks:
            result_df = pd.concat(chunks, ignore_index=True)
            
            # 如果需要，按每个RBP取top_n
            if top_n and 'rbp' in result_df.columns:
                result_df = get_top_n_per_rbp(result_df, top_n)
            
            # 保存结果
            result_df.to_csv(output_file, index=False)
        else:
            print("警告: 没有找到任何binding为True的记录")
            # 创建一个空的DataFrame保存
            result_df = pd.DataFrame(columns=sample_df.columns)
            result_df.to_csv(output_file, index=False)
    
    else:
        # 一次性读取小文件
        try:
            df = pd.read_csv(input_file)
            total_rows = len(df)
            
            # 过滤binding为True且分数大于阈值
            mask = df['binding'] == True
            if min_score > 0:
                mask = mask & (df['score'] >= min_score)
            
            result_df = df[mask].copy()
            filtered_rows = len(result_df)
            
            # 如果需要，按每个RBP取top_n
            if top_n and 'rbp' in result_df.columns:
                result_df = get_top_n_per_rbp(result_df, top_n)
            
            # 保存结果
            result_df.to_csv(output_file, index=False)
            
        except MemoryError:
            print("内存不足！尝试使用分块处理...")
            return filter_binding_true(input_file, output_file, min_score, top_n, chunk_size, verbose)
    
    if verbose:
        print(f"\n处理完成:")
        print(f"  总行数: {total_rows:,}")
        print(f"  保留行数: {filtered_rows:,}")
        print(f"  过滤比例: {filtered_rows/total_rows*100:.2f}%")
        print(f"  输出文件: {output_file}")
        
        # 显示一些统计信息
        if filtered_rows > 0 and 'rbp' in result_df.columns:
            print(f"\n按RBP统计 (前10个):")
            rbp_counts = result_df['rbp'].value_counts().head(10)
            for rbp, count in rbp_counts.items():
                print(f"  {rbp:<15}: {count:,}")
    
    return output_file, filtered_rows

def get_top_n_per_rbp(df, top_n):
    """
    对每个RBP取分数最高的top_n条记录
    """
    if 'score' not in df.columns or 'rbp' not in df.columns:
        return df
    
    # 按RBP分组，取每个RBP分数最高的top_n条记录
    result = df.groupby('rbp', group_keys=False).apply(
        lambda x: x.nlargest(min(top_n, len(x)), 'score')
    )
    
    return result

def filter_directory(input_dir, output_dir=None, 
                     min_score=0.0, top_n=None,
                     pattern="*.csv", recursive=False):
    """
    批量处理目录中的所有CSV文件
    """
    import glob
    
    input_path = Path(input_dir)
    
    # 确定输出目录
    if output_dir is None:
        output_dir = str(input_path.parent / f"{input_path.name}_filtered")
    os.makedirs(output_dir, exist_ok=True)
    
    # 查找CSV文件
    if recursive:
        file_pattern = f"**/{pattern}"
    else:
        file_pattern = pattern
    
    csv_files = list(input_path.glob(file_pattern))
    
    if not csv_files:
        print(f"在目录 {input_dir} 中没有找到匹配 {pattern} 的文件")
        return []
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    
    results = []
    for input_file in csv_files:
        # 跳过已经过滤过的文件
        if "_binding_true" in str(input_file) or "_filtered" in str(input_file):
            print(f"跳过可能已过滤的文件: {input_file.name}")
            continue
        
        output_file = os.path.join(output_dir, f"{input_file.stem}_binding_true.csv")
        
        print(f"\n处理: {input_file.name}")
        try:
            output_path, filtered_count = filter_binding_true(
                str(input_file), output_file, 
                min_score=min_score, top_n=top_n,
                verbose=False
            )
            
            if output_path:
                file_size_mb = os.path.getsize(output_path) / (1024 * 1024) if os.path.exists(output_path) else 0
                print(f"  保留 {filtered_count:,} 行，输出大小: {file_size_mb:.1f} MB")
                results.append((str(input_file), output_path, filtered_count))
        except Exception as e:
            print(f"  处理失败: {e}")
    
    return results

def analyze_filtered_results(input_file, detailed=False):
    """
    分析过滤后的结果文件
    """
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"无法读取文件: {e}")
        return
    
    print(f"文件: {input_file}")
    print(f"总行数: {len(df):,}")
    
    if len(df) == 0:
        print("文件为空")
        return
    
    # 基本统计
    if 'score' in df.columns:
        print(f"分数统计:")
        print(f"  平均值: {df['score'].mean():.4f}")
        print(f"  中位数: {df['score'].median():.4f}")
        print(f"  最大值: {df['score'].max():.4f}")
        print(f"  最小值: {df['score'].min():.4f}")
    
    # RBP统计
    if 'rbp' in df.columns:
        rbp_stats = df['rbp'].value_counts()
        print(f"\nRBP统计 (共 {len(rbp_stats)} 个不同RBP):")
        print(f"  Top 10 RBPs:")
        for rbp, count in rbp_stats.head(10).items():
            print(f"    {rbp:<15}: {count:,} ({count/len(df)*100:.1f}%)")
    
    # 模型类型统计
    if 'model_type' in df.columns:
        model_stats = df['model_type'].value_counts()
        print(f"\n模型类型统计:")
        for model_type, count in model_stats.items():
            print(f"  {model_type}: {count:,} ({count/len(df)*100:.1f}%)")
    
    # 详细统计
    if detailed and 'rbp' in df.columns and 'score' in df.columns:
        print(f"\n详细RBP统计:")
        for rbp in df['rbp'].unique()[:20]:  # 限制显示前20个
            rbp_data = df[df['rbp'] == rbp]
            if len(rbp_data) > 0:
                print(f"  {rbp:<15}: {len(rbp_data):>5,} 条, 平均分数: {rbp_data['score'].mean():.3f}")

def create_summary_report(filtered_files, output_file="filter_summary.csv"):
    """
    创建过滤结果摘要报告
    """
    if not filtered_files:
        print("没有可汇总的文件")
        return
    
    summary_data = []
    for input_file, output_file, filtered_count in filtered_files:
        # 获取原始文件大小
        original_size = os.path.getsize(input_file) if os.path.exists(input_file) else 0
        filtered_size = os.path.getsize(output_file) if os.path.exists(output_file) else 0
        
        # 尝试获取原始行数
        try:
            with open(input_file, 'r') as f:
                original_lines = sum(1 for _ in f) - 1  # 减去标题行
        except:
            original_lines = 0
        
        summary_data.append({
            'original_file': os.path.basename(input_file),
            'filtered_file': os.path.basename(output_file),
            'original_size_mb': original_size / (1024*1024),
            'filtered_size_mb': filtered_size / (1024*1024),
            'original_lines': original_lines,
            'filtered_lines': filtered_count,
            'reduction_ratio': (original_lines - filtered_count) / original_lines * 100 if original_lines > 0 else 0
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False, float_format='%.2f')
    
    print(f"\n过滤摘要已保存到: {output_file}")
    print(summary_df[['original_file', 'original_lines', 'filtered_lines', 'reduction_ratio']].to_string(index=False))

def main():
    parser = argparse.ArgumentParser(description='过滤RBP分析结果，保留binding为True的记录')
    parser.add_argument('input', help='输入CSV文件或包含CSV文件的目录')
    parser.add_argument('--output', '-o', help='输出文件或目录路径')
    parser.add_argument('--min-score', type=float, default=0.0, 
                       help='最小分数阈值 (默认: 0.0)')
    parser.add_argument('--top-n', type=int, 
                       help='每个RBP只保留分数最高的前N条记录')
    parser.add_argument('--batch', action='store_true',
                       help='批量处理目录中的所有CSV文件')
    parser.add_argument('--recursive', '-r', action='store_true',
                       help='递归查找子目录中的CSV文件（与--batch一起使用）')
    parser.add_argument('--pattern', default="*.csv",
                       help='文件匹配模式（默认: *.csv）')
    parser.add_argument('--analyze', action='store_true',
                       help='分析过滤后的结果文件')
    parser.add_argument('--detailed', action='store_true',
                       help='显示详细分析结果（与--analyze一起使用）')
    parser.add_argument('--chunk-size', type=int, default=100000,
                       help='分块处理大小（默认: 100000）')
    parser.add_argument('--summary', action='store_true',
                       help='创建过滤结果摘要报告')
    
    args = parser.parse_args()
    
    # 检查输入路径
    if not os.path.exists(args.input):
        print(f"错误: 输入路径不存在: {args.input}")
        return 1
    
    # 分析模式
    if args.analyze:
        analyze_filtered_results(args.input, args.detailed)
        return 0
    
    # 批量处理目录
    if args.batch or os.path.isdir(args.input):
        print("批量处理模式...")
        results = filter_directory(
            args.input, args.output,
            min_score=args.min_score,
            top_n=args.top_n,
            pattern=args.pattern,
            recursive=args.recursive
        )
        
        if args.summary and results:
            create_summary_report(results)
        
        return 0
    
    # 处理单个文件
    output_file, filtered_count = filter_binding_true(
        args.input, args.output,
        min_score=args.min_score,
        top_n=args.top_n,
        chunk_size=args.chunk_size,
        verbose=True
    )
    
    if output_file and args.summary:
        # 创建单个文件的摘要
        summary_data = [(args.input, output_file, filtered_count)]
        create_summary_report(summary_data, "filter_summary_single.csv")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
