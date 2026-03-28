#!/usr/bin/env python3
"""
AS分析主控脚本
用法：
    python run_as_pipeline.py --base_dir /path/to/base --output_dir /path/to/as
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from collections import defaultdict

# ========== 配置 ==========
CONDA_ENV = "common_env"
MIN_SAMPLES = 3          # 至少出现在3个样本中
MIN_TPM = 0.1            # 平均TPM > 0.1
# SUPPA 支持的事件类型：SE, SS (A5/A3), MX, RI, FL (AF/AL)
EVENT_TYPES = ['SE', 'SS', 'MX', 'RI', 'FL']

# ========== 函数定义 ==========

def run_command(cmd, conda=True):
    """运行命令，可选择是否在conda环境中运行"""
    if conda:
        full_cmd = ["conda", "run", "-n", CONDA_ENV] + cmd
    else:
        full_cmd = cmd
    print(f"Running: {' '.join(full_cmd)}")
    result = subprocess.run(full_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(full_cmd)}")
    return result.stdout

def build_tpm_matrix(base_dir, output_file):
    """从gffcompare生成的tmap文件中构建参考转录本TPM矩阵（SUPPA兼容格式）"""
    gffcompare_dir = os.path.join(base_dir, "gtf_merge", "gffcompare")
    if not os.path.isdir(gffcompare_dir):
        raise FileNotFoundError(f"gffcompare directory not found: {gffcompare_dir}")

    # 收集所有tmap文件
    tmap_files = []
    for sample_dir in os.listdir(gffcompare_dir):
        sample_path = os.path.join(gffcompare_dir, sample_dir)
        if os.path.isdir(sample_path):
            for f in os.listdir(sample_path):
                if f.endswith('.tmap') and 'ballgown' in f:
                    tmap_files.append((sample_dir, os.path.join(sample_path, f)))
                    break
    if not tmap_files:
        raise RuntimeError("No tmap files found in gffcompare directory")

    sample_data = defaultdict(lambda: defaultdict(float))
    all_ref_transcripts = set()

    for sample_name, tmap_file in tmap_files:
        print(f"Reading {tmap_file}...")
        try:
            df = pd.read_csv(tmap_file, sep='\t')
            # 筛选有效的参考转录本（ref_id != '-' 且 TPM > 0）
            valid = df[(df['ref_id'] != '-') & (df['TPM'] > 0)]
            for _, row in valid.iterrows():
                ref_id = row['ref_id']
                tpm = row['TPM']
                sample_data[sample_name][ref_id] += tpm
                all_ref_transcripts.add(ref_id)
        except Exception as e:
            print(f"Error reading {tmap_file}: {e}")

    # 构建DataFrame
    matrix = pd.DataFrame(index=list(all_ref_transcripts))
    for sample, trans_dict in sample_data.items():
        matrix[sample] = pd.Series(trans_dict)
    matrix = matrix.fillna(0)

    # 保存为SUPPA兼容格式：第一行是样本名（无空列），第一列是转录本ID（无标题）
    with open(output_file, 'w') as f:
        f.write('\t'.join(matrix.columns) + '\n')
    matrix.to_csv(output_file, sep='\t', header=False, index=True, mode='a')
    print(f"Transcript TPM matrix saved to {output_file}, shape={matrix.shape}")
    return matrix

def filter_expressed_gtf(ref_gtf, tpm_matrix, expressed_gtf):
    """从参考注释中筛选表达的转录本，生成expressed GTF"""
    # 读取TPM矩阵，确定表达的转录本（在至少MIN_SAMPLES个样本中TPM > MIN_TPM）
    tpm_df = pd.read_csv(tpm_matrix, sep='\t', index_col=0)
    expressed_in_samples = (tpm_df > MIN_TPM).sum(axis=1)
    expressed_transcripts = expressed_in_samples[expressed_in_samples >= MIN_SAMPLES].index.tolist()
    print(f"Expressed transcripts: {len(expressed_transcripts)}")

    # 从参考GTF中提取这些转录本的行
    output_lines = []
    with open(ref_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                output_lines.append(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                output_lines.append(line)
                continue
            # 提取transcript_id
            attrs = parts[8]
            tid = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tid = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break
            if tid and tid in expressed_transcripts:
                output_lines.append(line)

    with open(expressed_gtf, 'w') as f:
        f.writelines(output_lines)
    print(f"Expressed GTF written to {expressed_gtf}")

def extract_high_quality_novel(base_dir, output_file):
    """从tmap文件中提取高质量新转录本（class_code='u'），输出详细TSV"""
    gffcompare_dir = os.path.join(base_dir, "gtf_merge", "gffcompare")
    if not os.path.isdir(gffcompare_dir):
        raise FileNotFoundError(f"gffcompare directory not found: {gffcompare_dir}")

    transcript_counts = defaultdict(int)
    transcript_info = defaultdict(dict)  # 存储每个样本中的信息

    for sample_dir in os.listdir(gffcompare_dir):
        sample_path = os.path.join(gffcompare_dir, sample_dir)
        if not os.path.isdir(sample_path):
            continue
        # 查找tmap文件
        tmap_file = None
        for f in os.listdir(sample_path):
            if f.endswith('.tmap') and 'ballgown' in f:
                tmap_file = os.path.join(sample_path, f)
                break
        if not tmap_file:
            continue

        print(f"Reading {tmap_file}...")
        try:
            df = pd.read_csv(tmap_file, sep='\t')
            novel_df = df[df['class_code'] == 'u']
            for _, row in novel_df.iterrows():
                qry_id = row['qry_id']
                if qry_id != '-':
                    transcript_counts[qry_id] += 1
                    transcript_info[qry_id][sample_dir] = {
                        'TPM': row['TPM'],
                        'cov': row['cov'],
                        'len': row['len']
                    }
        except Exception as e:
            print(f"Error reading {tmap_file}: {e}")

    # 筛选高质量新转录本
    output_data = []
    for tid, count in transcript_counts.items():
        if count >= MIN_SAMPLES:
            samples_info = transcript_info[tid]
            tpm_values = [info['TPM'] for info in samples_info.values()]
            cov_values = [info['cov'] for info in samples_info.values()]
            len_values = [info['len'] for info in samples_info.values()]
            avg_tpm = np.mean(tpm_values)
            if avg_tpm > MIN_TPM:
                output_data.append({
                    'transcript_id': tid,
                    'sample_count': count,
                    'avg_TPM': avg_tpm,
                    'max_TPM': max(tpm_values),
                    'avg_coverage': np.mean(cov_values),
                    'avg_length': np.mean(len_values),
                    'samples': ','.join(samples_info.keys())
                })

    # 转换为DataFrame并保存
    df_output = pd.DataFrame(output_data)
    df_output = df_output.sort_values(['sample_count', 'avg_TPM'], ascending=[False, False])
    df_output.to_csv(output_file, sep='\t', index=False)
    print(f"High-quality novel transcripts saved to {output_file}")
    print(f"  Total: {len(df_output)} transcripts")
    return df_output

def create_enhanced_annotation(expressed_gtf, merged_gtf, novel_tsv, enhanced_gtf):
    """将高质量新转录本从merged.gtf追加到expressed.gtf，生成enhanced GTF"""
    # 读取高质量转录本ID（从TSV的第一列）
    if not os.path.exists(novel_tsv):
        print(f"Warning: {novel_tsv} not found, no novel transcripts to add.")
        high_quality_ids = set()
    else:
        df_novel = pd.read_csv(novel_tsv, sep='\t')
        high_quality_ids = set(df_novel['transcript_id'].tolist())
    print(f"Adding {len(high_quality_ids)} transcripts to enhanced GTF")

    # 复制expressed GTF到输出
    with open(expressed_gtf, 'r') as f_in, open(enhanced_gtf, 'w') as f_out:
        f_out.writelines(f_in.readlines())

    # 从merged.gtf中提取高质量转录本的条目
    with open(merged_gtf, 'r') as f_in:
        current_tid = None
        buffer = []
        for line in f_in:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            attrs = parts[8]
            tid = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tid = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break
            if tid:
                if tid != current_tid:
                    # 上一个转录本结束，检查是否保存
                    if current_tid and current_tid in high_quality_ids:
                        with open(enhanced_gtf, 'a') as f_out:
                            f_out.writelines(buffer)
                    current_tid = tid
                    buffer = [line]
                else:
                    buffer.append(line)
        # 处理最后一个
        if current_tid and current_tid in high_quality_ids:
            with open(enhanced_gtf, 'a') as f_out:
                f_out.writelines(buffer)

    print(f"Enhanced GTF written to {enhanced_gtf}")

def generate_events(gtf_file, output_dir):
    """使用SUPPA生成事件文件"""
    os.makedirs(output_dir, exist_ok=True)
    for event_type in EVENT_TYPES:
        out_prefix = os.path.join(output_dir, f"{event_type}")
        cmd = ["suppa.py", "generateEvents", "-i", gtf_file, "-o", out_prefix, "-e", event_type, "-f", "ioe"]
        run_command(cmd, conda=True)
    print(f"Events generated in {output_dir}")

def compute_psi(events_dir, tpm_matrix, output_dir):
    """使用SUPPA计算PSI值"""
    os.makedirs(output_dir, exist_ok=True)
    # 获取所有事件文件（.ioe）
    ioe_files = [f for f in os.listdir(events_dir) if f.endswith('.ioe')]
    for ioe in ioe_files:
        # 提取事件类型（去掉_strict.ioe后缀，但实际文件可能有_strict）
        event_base = ioe.replace('.ioe', '')
        event_file = os.path.join(events_dir, ioe)
        out_prefix = os.path.join(output_dir, event_base)
        cmd = ["suppa.py", "psiPerEvent", "-i", event_file, "-e", tpm_matrix, "-o", out_prefix]
        run_command(cmd, conda=True)
    print(f"PSI values computed in {output_dir}")

# ========== 主函数 ==========
def main():
    parser = argparse.ArgumentParser(description="AS analysis pipeline")
    parser.add_argument("--base_dir", required=True, help="Base working directory (where gtf_merge, etc. are located)")
    parser.add_argument("--output_dir", required=True, help="Output directory for AS results (e.g., ./as)")
    args = parser.parse_args()

    # 确定路径
    base_dir = os.path.abspath(args.base_dir)
    output_dir = os.path.abspath(args.output_dir)
    merged_gtf = os.path.join(base_dir, "gtf_merge", "merged.gtf")
    ref_gtf = os.path.join(base_dir, "gencode.v45.annotation.gtf")

    # 创建输出目录
    expressed_dir = os.path.join(output_dir, "expressed")
    enhanced_dir = os.path.join(output_dir, "enhanced")
    for d in [expressed_dir, enhanced_dir]:
        os.makedirs(d, exist_ok=True)

    # 中间文件路径
    tpm_matrix = os.path.join(output_dir, "transcript_tpm_matrix.tsv")
    expressed_gtf = os.path.join(output_dir, "gencode.v45.expressed.gtf")
    enhanced_gtf = os.path.join(output_dir, "gencode.v45.enhanced.gtf")
    novel_tsv = os.path.join(output_dir, "high_quality_novel_transcripts.tsv")

    # Step 1: 构建TPM矩阵（从tmap文件）
    print("="*50)
    print("Step 1: Building TPM matrix from tmap files")
    tpm_df = build_tpm_matrix(base_dir, tpm_matrix)

    # Step 2: 创建expressed GTF
    print("="*50)
    print("Step 2: Creating expressed GTF")
    filter_expressed_gtf(ref_gtf, tpm_matrix, expressed_gtf)

    # Step 3: 提取高质量新转录本（输出详细TSV）
    print("="*50)
    print("Step 3: Extracting high-quality novel transcripts")
    extract_high_quality_novel(base_dir, novel_tsv)

    # Step 4: 创建enhanced GTF
    print("="*50)
    print("Step 4: Creating enhanced GTF")
    create_enhanced_annotation(expressed_gtf, merged_gtf, novel_tsv, enhanced_gtf)

    # Step 5: 为两组注释生成事件文件
    print("="*50)
    print("Step 5: Generating events for expressed group")
    generate_events(expressed_gtf, os.path.join(expressed_dir, "events"))

    print("="*50)
    print("Step 5b: Generating events for enhanced group")
    generate_events(enhanced_gtf, os.path.join(enhanced_dir, "events"))

    # Step 6: 计算PSI
    print("="*50)
    print("Step 6: Computing PSI for expressed group")
    compute_psi(os.path.join(expressed_dir, "events"), tpm_matrix, os.path.join(expressed_dir, "psi"))

    print("="*50)
    print("Step 6b: Computing PSI for enhanced group")
    compute_psi(os.path.join(enhanced_dir, "events"), tpm_matrix, os.path.join(enhanced_dir, "psi"))

    print("="*50)
    print("AS analysis completed successfully!")

if __name__ == "__main__":
    main()