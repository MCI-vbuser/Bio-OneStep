#!/usr/bin/env python
"""
高级RBP分析功能
"""
import pandas as pd
import numpy as np
from rbp_inference import RBPInference
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns

class AdvancedRBPAnalysis:
    def __init__(self, model_dir="../Results/PARCLIP_models/"):
        self.inferencer = RBPInference(model_dir)
    
    def analyze_co_binding_patterns(self, results_df, threshold=0.7):
        """
        分析RBP共结合模式
        """
        # 筛选高置信度结合
        high_conf = results_df[results_df['score'] > threshold]
        
        # 为每个序列记录结合的RBPs
        sequence_rbps = defaultdict(set)
        for _, row in high_conf.iterrows():
            if row['binding']:
                sequence_rbps[row['sequence']].add(row['rbp'])
        
        # 计算RBP共出现频率
        co_occurrence = defaultdict(int)
        for rbps in sequence_rbps.values():
            rbp_list = list(rbps)
            for i in range(len(rbp_list)):
                for j in range(i+1, len(rbp_list)):
                    pair = tuple(sorted([rbp_list[i], rbp_list[j]]))
                    co_occurrence[pair] += 1
        
        # 转换为DataFrame
        co_df = pd.DataFrame([
            {'rbp1': p[0], 'rbp2': p[1], 'co_occurrence': count}
            for p, count in co_occurrence.items()
        ])
        
        if not co_df.empty:
            co_df = co_df.sort_values('co_occurrence', ascending=False)
        
        return co_df, dict(sequence_rbps)
    
    def build_rbp_network(self, results_df, min_co_occurrence=3):
        """
        构建RBP相互作用网络
        """
        # 分析共结合模式
        co_df, _ = self.analyze_co_binding_patterns(results_df)
        
        # 创建网络图
        G = nx.Graph()
        
        # 添加节点（所有RBPs）
        all_rbps = results_df['rbp'].unique()
        for rbp in all_rbps:
            G.add_node(rbp)
        
        # 添加边（共结合关系）
        if not co_df.empty:
            for _, row in co_df[co_df['co_occurrence'] >= min_co_occurrence].iterrows():
                G.add_edge(row['rbp1'], row['rbp2'], 
                          weight=row['co_occurrence'])
        
        return G
    
    def plot_rbp_network(self, results_df, save_path=None, min_co_occurrence=3):
        """
        绘制RBP相互作用网络图
        """
        # 构建网络
        G = self.build_rbp_network(results_df, min_co_occurrence)
        
        if len(G.edges()) == 0:
            print("No co-binding relationships found above threshold.")
            return None
        
        # 计算节点属性
        node_sizes = []
        node_colors = []
        
        # 计算每个RBP的结合频率作为节点大小
        binding_freq = results_df.groupby('rbp')['binding'].mean()
        
        for node in G.nodes():
            freq = binding_freq.get(node, 0)
            node_sizes.append(1000 + freq * 5000)  # 根据频率调整大小
            node_colors.append(freq)  # 用频率作为颜色
        
        # 创建图表
        plt.figure(figsize=(16, 12))
        
        # 使用spring布局
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # 绘制边
        edges = G.edges(data=True)
        edge_weights = [data['weight'] for _, _, data in edges]
        edge_widths = [w / max(edge_weights) * 5 for w in edge_weights]
        
        nx.draw_networkx_edges(G, pos, 
                              width=edge_widths,
                              alpha=0.6,
                              edge_color='gray')
        
        # 绘制节点
        nodes = nx.draw_networkx_nodes(G, pos,
                                      node_size=node_sizes,
                                      node_color=node_colors,
                                      cmap=plt.cm.YlOrRd,
                                      alpha=0.8)
        
        # 绘制标签
        nx.draw_networkx_labels(G, pos, 
                               font_size=10,
                               font_weight='bold')
        
        # 添加颜色条
        sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd, 
                                  norm=plt.Normalize(vmin=min(node_colors), 
                                                    vmax=max(node_colors)))
        sm.set_array([])
        plt.colorbar(sm, label='Binding Frequency')
        
        plt.title('RBP Co-binding Network\n(Node size and color indicate binding frequency)',
                 fontsize=14, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Network plot saved to: {save_path}")
        
        plt.show()
        
        # 输出网络统计信息
        print(f"Network Statistics:")
        print(f"  Nodes: {G.number_of_nodes()}")
        print(f"  Edges: {G.number_of_edges()}")
        print(f"  Average degree: {np.mean([d for n, d in G.degree()]):.2f}")
        print(f"  Network density: {nx.density(G):.3f}")
        
        # 计算中心性指标
        if G.number_of_nodes() > 0:
            degree_centrality = nx.degree_centrality(G)
            betweenness_centrality = nx.betweenness_centrality(G)
            
            print("\nTop 5 RBPs by Degree Centrality:")
            for rbp, centrality in sorted(degree_centrality.items(), 
                                        key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {rbp}: {centrality:.3f}")
            
            print("\nTop 5 RBPs by Betweenness Centrality:")
            for rbp, centrality in sorted(betweenness_centrality.items(), 
                                        key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {rbp}: {centrality:.3f}")
        
        return G
    
    def analyze_sequence_features(self, sequences, rbp_scores):
        """
        分析序列特征与RBP结合的关系
        """
        from Bio.SeqUtils import GC
        
        features = []
        
        for seq, scores in zip(sequences, rbp_scores):
            # 计算序列特征
            gc_content = GC(seq)
            length = len(seq)
            
            # 计算平均结合分数
            avg_score = np.mean(scores) if len(scores) > 0 else 0
            max_score = np.max(scores) if len(scores) > 0 else 0
            
            features.append({
                'sequence': seq[:50] + '...' if len(seq) > 50 else seq,
                'length': length,
                'gc_content': gc_content,
                'avg_rbp_score': avg_score,
                'max_rbp_score': max_score,
                'num_binding_rbps': sum(1 for s in scores if s > 0.5)
            })
        
        return pd.DataFrame(features)
    
    def correlate_with_expression(self, rbp_results, expression_data):
        """
        将RBP结合与基因表达数据进行关联分析
        """
        # 这里需要根据具体数据格式进行调整
        # 假设expression_data包含gene_name和expression_value
        
        correlation_results = []
        
        for rbp in rbp_results['rbp'].unique():
            rbp_data = rbp_results[rbp_results['rbp'] == rbp]
            
            # 合并表达数据
            merged = pd.merge(rbp_data, expression_data,
                            left_on='sequence', right_on='gene_name',
                            how='inner')
            
            if len(merged) > 10:  # 需要有足够的数据点
                # 计算相关性
                correlation = merged['score'].corr(merged['expression_value'])
                
                correlation_results.append({
                    'rbp': rbp,
                    'correlation': correlation,
                    'n_samples': len(merged),
                    'avg_score': merged['score'].mean(),
                    'avg_expression': merged['expression_value'].mean()
                })
        
        return pd.DataFrame(correlation_results).sort_values('correlation', 
                                                           key=abs, 
                                                           ascending=False)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='高级RBP分析')
    parser.add_argument('--results', type=str, required=True,
                       help='结果CSV文件路径')
    parser.add_argument('--output', type=str, default='./advanced_analysis',
                       help='输出目录')
    parser.add_argument('--network', action='store_true',
                       help='构建和可视化RBP网络')
    parser.add_argument('--co-binding', action='store_true',
                       help='分析共结合模式')
    
    args = parser.parse_args()
    
    # 读取结果
    results_df = pd.read_csv(args.results)
    
    # 初始化分析器
    analyzer = AdvancedRBPAnalysis()
    
    # 创建输出目录
    import os
    os.makedirs(args.output, exist_ok=True)
    
    if args.co_binding:
        print("Analyzing co-binding patterns...")
        co_df, seq_rbps = analyzer.analyze_co_binding_patterns(results_df)
        
        co_output = os.path.join(args.output, 'co_binding_analysis.csv')
        co_df.to_csv(co_output, index=False)
        print(f"Co-binding analysis saved to: {co_output}")
        
        # 显示Top共结合对
        print("\nTop 10 co-binding RBP pairs:")
        for _, row in co_df.head(10).iterrows():
            print(f"  {row['rbp1']} - {row['rbp2']}: {row['co_occurrence']} sequences")
    
    if args.network:
        print("\nBuilding RBP interaction network...")
        network_output = os.path.join(args.output, 'rbp_network.png')
        analyzer.plot_rbp_network(results_df, save_path=network_output)
