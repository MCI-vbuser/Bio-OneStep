#!/usr/bin/env python
"""
DeepRiPe RBP结合预测推理脚本
支持批量处理自定义序列
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 添加自定义模块路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from deeptripe_utils import *
from integrated_gradients import IntegratedGradients
from tensorflow.keras.models import load_model

# RBP名称定义（与原模型一致）
RBPNAMES = {
    'high': np.array(['DND1', 'CPSF7', 'CPSF6', 'CPSF1', 'CSTF2', 'CSTF2T', 
                     'ZC3H7B', 'FMR1iso1', 'RBM10', 'MOV10', 'ELAVL1']),
    'med': np.array(['TARDBP', 'ELAVL2', 'ELAVL3', 'ELAVL4', 'RBM20', 
                    'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'EWSR1', 'HNRNPD', 
                    'RBPMS', 'SRRM4', 'AGO2', 'NUDT21', 'FIP1L1', 'CAPRIN1', 
                    'FMR1iso7', 'FXR2', 'AGO1', 'L1RE1', 'ORF1']),
    'low': np.array(['MBNL1', 'P53_NONO', 'PUM2', 'QKI', 'AGO3', 'FUS', 
                    'TAF15', 'ZFP36', 'DICER1', 'EIF3A', 'EIF3D', 'EIF3G', 
                    'SSB', 'PAPD5', 'CPSF4', 'CPSF3', 'RTCB', 'FXR1', 
                    'NOP58', 'NOP56', 'FBL', 'LIN28A', 'LIN28B', 'UPF1', 
                    'G35', 'G45', 'XPO5'])
}

class RBPInference:
    def __init__(self, model_dir="../Results/PARCLIP_models/"):
        """
        初始化RBP推理器
        """
        self.model_dir = Path(model_dir)
        self.models = {}
        self.ig_res = {}
        
        # 加载所有模型
        self._load_models()
    
    def _load_models(self):
        """加载预训练模型"""
        model_files = {
            'high': self.model_dir / 'model_RBPshigh.h5',
            'med': self.model_dir / 'model_RBPsmed.h5',
            'low': self.model_dir / 'model_RBPslow.h5'
        }
        
        custom_objects = {'precision': precision, 'recall': recall}
        
        for model_type, model_path in model_files.items():
            if model_path.exists():
                print(f"Loading {model_type} model...")
                self.models[model_type] = load_model(
                    model_path, 
                    custom_objects=custom_objects
                )
                self.ig_res[model_type] = IntegratedGradients(self.models[model_type])
            else:
                print(f"Warning: Model file not found: {model_path}")
    
    def predict_sequences(self, sequences, sequence_names=None):
        """
        预测RNA序列的RBP结合
        
        Parameters:
        -----------
        sequences : list of str
            RNA序列列表
        sequence_names : list of str, optional
            序列名称列表
        
        Returns:
        --------
        pandas.DataFrame
            预测结果DataFrame
        """
        if sequence_names is None:
            sequence_names = [f"seq_{i}" for i in range(len(sequences))]
        
        # 准备输入数据
        X_seq = np.array([seq_to_mat(seq) for seq in sequences])
        X_region = np.full((len(sequences), 250, 4), 0.25)  # 均匀区域输入
        
        results = []
        
        # 对每个模型进行预测
        for model_type in ['high', 'med', 'low']:
            if model_type in self.models:
                model = self.models[model_type]
                predictions = model.predict([X_seq, X_region])
                
                # 为每个序列和每个RBP保存结果
                for seq_idx, seq_name in enumerate(sequence_names):
                    for rbp_idx, rbp_name in enumerate(RBPNAMES[model_type]):
                        score = predictions[seq_idx, rbp_idx]
                        results.append({
                            'sequence': seq_name,
                            'rbp': rbp_name,
                            'model_type': model_type,
                            'score': float(score),
                            'threshold': 0.5,  # 可根据需要调整
                            'binding': score > 0.5
                        })
        
        return pd.DataFrame(results)
    
    def predict_bed_regions(self, bed_file, genome_fasta, flank=75):
        """
        从BED文件预测RBP结合
        
        Parameters:
        -----------
        bed_file : str
            BED文件路径
        genome_fasta : str
            基因组FASTA文件路径
        flank : int
            上下游延伸长度
        
        Returns:
        --------
        pandas.DataFrame
            预测结果
        """
        try:
            import pybedtools
        except ImportError:
            print("Please install pybedtools: pip install pybedtools")
            return None
        
        # 读取BED文件
        bed = pybedtools.BedTool(bed_file)
        
        # 提取序列
        sequences = []
        seq_names = []
        
        for interval in bed:
            # 扩展区域
            start = max(0, interval.start - flank)
            end = interval.end + flank
            
            # 创建临时bed对象
            temp_bed = pybedtools.BedTool(
                f"{interval.chrom}\t{start}\t{end}\t{interval.name}\t.\t{interval.strand}",
                from_string=True
            )
            
            # 提取序列
            temp_fasta = temp_bed.sequence(fi=genome_fasta, s=True)
            
            with open(temp_fasta.seqfn) as f:
                for line in f:
                    if line.startswith('>'):
                        seq_names.append(line.strip()[1:])
                    else:
                        sequences.append(line.strip().upper().replace('T', 'U'))
        
        # 进行预测
        return self.predict_sequences(sequences, seq_names)
    
    def get_top_rbps(self, predictions_df, top_n=10):
        """
        获取结合分数最高的RBP
        
        Parameters:
        -----------
        predictions_df : pandas.DataFrame
            预测结果DataFrame
        top_n : int
            返回前N个结果
        
        Returns:
        --------
        pandas.DataFrame
            排序后的结果
        """
        return (predictions_df
                .sort_values('score', ascending=False)
                .head(top_n))
    
    def analyze_novel_transcripts(self, fasta_file, output_csv=None):
        """
        分析新转录本的RBP结合模式
        
        Parameters:
        -----------
        fasta_file : str
            新转录本FASTA文件路径
        output_csv : str, optional
            输出CSV文件路径
        
        Returns:
        --------
        pandas.DataFrame
            分析结果
        """
        from Bio import SeqIO
        
        sequences = []
        seq_names = []
        
        # 读取FASTA文件
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_names.append(record.id)
            sequences.append(str(record.seq).upper().replace('T', 'U'))
        
        # 进行预测
        results = self.predict_sequences(sequences, seq_names)
        
        if output_csv:
            results.to_csv(output_csv, index=False)
            print(f"Results saved to {output_csv}")
        
        return results

def main():
    parser = argparse.ArgumentParser(description='DeepRiPe RBP结合预测')
    parser.add_argument('--input', type=str, required=True,
                       help='输入文件：FASTA或BED格式')
    parser.add_argument('--input-type', choices=['fasta', 'bed'], required=True,
                       help='输入文件类型')
    parser.add_argument('--genome', type=str, 
                       help='基因组FASTA文件（BED输入时需要）')
    parser.add_argument('--output', type=str, default='rbp_predictions.csv',
                       help='输出CSV文件路径')
    parser.add_argument('--model-dir', type=str, default='../Results/PARCLIP_models/',
                       help='模型目录路径')
    parser.add_argument('--top-n', type=int, default=20,
                       help='显示前N个结果')
    
    args = parser.parse_args()
    
    # 初始化推理器
    print("Initializing RBP Inference Engine...")
    inferencer = RBPInference(model_dir=args.model_dir)
    
    # 根据输入类型进行处理
    if args.input_type == 'fasta':
        print(f"Analyzing FASTA file: {args.input}")
        results = inferencer.analyze_novel_transcripts(
            args.input, 
            args.output
        )
    elif args.input_type == 'bed':
        if not args.genome:
            parser.error("--genome required for BED input")
        print(f"Analyzing BED regions: {args.input}")
        results = inferencer.predict_bed_regions(
            args.input, 
            args.genome
        )
    
    # 保存结果
    if results is not None:
        results.to_csv(args.output, index=False)
        
        # 显示前N个结果
        print("\n" + "="*60)
        print(f"Top {args.top_n} RBP binding predictions:")
        print("="*60)
        
        top_results = inferencer.get_top_rbps(results, args.top_n)
        for idx, row in top_results.iterrows():
            print(f"{row['rbp']:15s} | {row['sequence']:30s} | "
                  f"Score: {row['score']:.3f} | Binding: {row['binding']}")
        
        print(f"\nComplete results saved to: {args.output}")
        
        # 统计信息
        binding_summary = results.groupby('rbp')['binding'].sum().sort_values(ascending=False)
        print("\nRBP binding frequency:")
        for rbp, count in binding_summary.head(10).items():
            total = len(results[results['rbp'] == rbp])
            print(f"  {rbp:15s}: {count}/{total} ({count/total:.1%})")

if __name__ == "__main__":
    main()
