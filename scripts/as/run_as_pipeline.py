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
MIN_SAMPLES = 3
MIN_TPM = 0.1
EVENT_TYPES = ['SE', 'SS', 'MX', 'RI', 'FL']

def run_command(cmd, conda=True):
    """运行命令，实时输出（对于 python 脚本强制 -u）"""
    if conda and cmd[0] == "python":
        # 如果命令是 python script.py ...，则插入 -u 确保无缓冲
        cmd = ["python", "-u"] + cmd[1:]
    if conda:
        full_cmd = ["conda", "run", "-n", CONDA_ENV] + cmd
    else:
        full_cmd = cmd
    print(f"Running: {' '.join(full_cmd)}")
    process = subprocess.Popen(full_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, bufsize=1)
    for line in process.stdout:
        print(line, end='')
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(full_cmd)}")

def build_tpm_matrix(base_dir, output_file):
    """从tmap文件中构建TPM矩阵，包含参考转录本和新转录本（class_code='u'）"""
    gffcompare_dir = os.path.join(base_dir, "gtf_merge", "gffcompare")
    if not os.path.isdir(gffcompare_dir):
        raise FileNotFoundError(f"gffcompare directory not found: {gffcompare_dir}")

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
    all_transcripts = set()  # 包括参考转录本和新转录本

    for sample_name, tmap_file in tmap_files:
        print(f"Reading {tmap_file}...")
        try:
            df = pd.read_csv(tmap_file, sep='\t')
            # 处理参考转录本 (ref_id != '-')
            ref_valid = df[(df['ref_id'] != '-') & (df['TPM'] > 0)]
            for _, row in ref_valid.iterrows():
                ref_id = row['ref_id']
                tpm = row['TPM']
                sample_data[sample_name][ref_id] += tpm
                all_transcripts.add(ref_id)

            # 处理新转录本 (class_code='u' 且 qry_id != '-')
            novel_valid = df[(df['class_code'] == 'u') & (df['qry_id'] != '-') & (df['TPM'] > 0)]
            for _, row in novel_valid.iterrows():
                qry_id = row['qry_id']
                tpm = row['TPM']
                sample_data[sample_name][qry_id] += tpm
                all_transcripts.add(qry_id)
        except Exception as e:
            print(f"Error reading {tmap_file}: {e}")

    # 构建DataFrame
    matrix = pd.DataFrame(index=list(all_transcripts))
    for sample, trans_dict in sample_data.items():
        matrix[sample] = pd.Series(trans_dict)
    matrix = matrix.fillna(0)

    # 保存为SUPPA兼容格式
    with open(output_file, 'w') as f:
        f.write('\t'.join(matrix.columns) + '\n')
    matrix.to_csv(output_file, sep='\t', header=False, index=True, mode='a')
    print(f"Transcript TPM matrix saved to {output_file}, shape={matrix.shape}")
    return matrix

def filter_expressed_gtf(ref_gtf, tpm_matrix, expressed_gtf):
    """从参考注释中筛选表达的转录本（基于TPM矩阵中的转录本ID）"""
    tpm_df = pd.read_csv(tpm_matrix, sep='\t', index_col=0)
    # 矩阵中的行索引即为所有转录本ID（包括参考和新转录本），但这里我们只保留参考注释中存在的转录本
    expressed_transcripts = set(tpm_df.index)
    print(f"Total transcripts in TPM matrix: {len(expressed_transcripts)}")

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
    print(f"Expressed GTF written to {expressed_gtf}, containing {len(expressed_transcripts)} transcripts")

def extract_high_quality_novel(base_dir, output_file):
    """提取高质量新转录本（与之前相同，无需修改）"""
    gffcompare_dir = os.path.join(base_dir, "gtf_merge", "gffcompare")
    if not os.path.isdir(gffcompare_dir):
        raise FileNotFoundError(f"gffcompare directory not found: {gffcompare_dir}")

    transcript_counts = defaultdict(int)
    transcript_info = defaultdict(dict)

    for sample_dir in os.listdir(gffcompare_dir):
        sample_path = os.path.join(gffcompare_dir, sample_dir)
        if not os.path.isdir(sample_path):
            continue
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

    df_output = pd.DataFrame(output_data)
    df_output = df_output.sort_values(['sample_count', 'avg_TPM'], ascending=[False, False])
    df_output.to_csv(output_file, sep='\t', index=False)
    print(f"High-quality novel transcripts saved to {output_file}, total: {len(df_output)}")
    return df_output

def create_enhanced_annotation(expressed_gtf, merged_gtf, novel_tsv, enhanced_gtf):
    """将高质量新转录本从merged.gtf追加到expressed.gtf，生成enhanced GTF"""
    if not os.path.exists(novel_tsv):
        print(f"Warning: {novel_tsv} not found, no novel transcripts to add.")
        high_quality_ids = set()
    else:
        df_novel = pd.read_csv(novel_tsv, sep='\t')
        high_quality_ids = set(df_novel['transcript_id'].tolist())
    print(f"Adding {len(high_quality_ids)} transcripts to enhanced GTF")

    # 复制expressed GTF
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
                    if current_tid and current_tid in high_quality_ids:
                        with open(enhanced_gtf, 'a') as f_out:
                            f_out.writelines(buffer)
                    current_tid = tid
                    buffer = [line]
                else:
                    buffer.append(line)
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

def compute_psi(events_dir, tpm_matrix, output_dir):
    """使用SUPPA计算PSI值"""
    os.makedirs(output_dir, exist_ok=True)
    ioe_files = [f for f in os.listdir(events_dir) if f.endswith('.ioe')]
    for ioe in ioe_files:
        event_base = ioe.replace('.ioe', '')
        event_file = os.path.join(events_dir, ioe)
        out_prefix = os.path.join(output_dir, event_base)
        cmd = ["suppa.py", "psiPerEvent", "-i", event_file, "-e", tpm_matrix, "-o", out_prefix]
        run_command(cmd, conda=True)

def main():
    parser = argparse.ArgumentParser(description="AS analysis pipeline")
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    output_dir = os.path.abspath(args.output_dir)
    merged_gtf = os.path.join(base_dir, "gtf_merge", "merged.gtf")
    ref_gtf = os.path.join(base_dir, "gencode.v45.annotation.gtf")

    expressed_dir = os.path.join(output_dir, "expressed")
    enhanced_dir = os.path.join(output_dir, "enhanced")
    for d in [expressed_dir, enhanced_dir]:
        os.makedirs(d, exist_ok=True)

    tpm_matrix = os.path.join(output_dir, "transcript_tpm_matrix.tsv")
    expressed_gtf = os.path.join(output_dir, "gencode.v45.expressed.gtf")
    enhanced_gtf = os.path.join(output_dir, "gencode.v45.enhanced.gtf")
    novel_tsv = os.path.join(output_dir, "high_quality_novel_transcripts.tsv")

    print("="*50)
    print("Step 1: Building TPM matrix from tmap files")
    tpm_df = build_tpm_matrix(base_dir, tpm_matrix)

    print("="*50)
    print("Step 2: Creating expressed GTF")
    filter_expressed_gtf(ref_gtf, tpm_matrix, expressed_gtf)

    print("="*50)
    print("Step 3: Extracting high-quality novel transcripts")
    extract_high_quality_novel(base_dir, novel_tsv)

    print("="*50)
    print("Step 4: Creating enhanced GTF")
    create_enhanced_annotation(expressed_gtf, merged_gtf, novel_tsv, enhanced_gtf)

    print("="*50)
    print("Step 5: Generating events for expressed group")
    generate_events(expressed_gtf, os.path.join(expressed_dir, "events"))

    print("="*50)
    print("Step 5b: Generating events for enhanced group")
    generate_events(enhanced_gtf, os.path.join(enhanced_dir, "events"))

    print("="*50)
    print("Step 6: Computing PSI for expressed group")
    compute_psi(os.path.join(expressed_dir, "events"), tpm_matrix, os.path.join(expressed_dir, "psi"))

    print("="*50)
    print("Step 6b: Computing PSI for enhanced group")
    compute_psi(os.path.join(enhanced_dir, "events"), tpm_matrix, os.path.join(enhanced_dir, "psi"))

    print("="*50)
    print("AS analysis completed successfully!")

    print("="*50)
    print("Step 7: Generating statistics and report")
    import as_stats

    stats_list = []
    for i, dir_path in enumerate([expressed_dir, enhanced_dir]):
        group_name = "expressed" if i == 0 else "enhanced"
        try:
            stats = as_stats.process_group(dir_path, group_name)
            stats_list.append(stats)
        except Exception as e:
            print(f"Error processing {dir_path}: {e}")
            sys.exit(1)

    # 只为 expressed 组保存详细数据，传入 TPM 矩阵路径
    as_stats.save_detailed_psi_to_db(stats_list[1], tpm_matrix, os.path.join(output_dir, "stats.db"))
    as_stats.generate_html(stats_list, os.path.join(output_dir, "index.html"))
    print("Statistics report generated.")

if __name__ == "__main__":
    main()