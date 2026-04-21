"""
DeepRiPe工具函数 - 清理版本
"""
import numpy as np
import h5py
from tensorflow.keras import backend as K
import sys
import os

def precision(y_true, y_pred):
    """自定义精确率指标"""
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def recall(y_true, y_pred):
    """自定义召回率指标"""
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def seq_to_mat(seq, length=150):
    """
    将RNA序列转换为one-hot编码矩阵
    A->[1,0,0,0], C->[0,1,0,0], G->[0,0,1,0], U->[0,0,0,1]
    """
    seq = seq.upper()
    vocab = {'A':0, 'C':1, 'G':2, 'U':3}
    
    # 处理序列长度
    if len(seq) < length:
        seq = seq + 'N' * (length - len(seq))
    elif len(seq) > length:
        seq = seq[:length]
    
    mat = np.zeros((length, 4))
    for i, base in enumerate(seq):
        if base in vocab:
            mat[i, vocab[base]] = 1
        else:
            # 对于未知碱基，使用均匀分布
            mat[i, :] = 0.25
    return mat

def seq_to_1hot(seq):
    """与seq_to_mat相同，但保持接口兼容性"""
    return seq_to_mat(seq)

def load_data(h5_path):
    """
    从h5文件加载数据
    """
    with h5py.File(h5_path, 'r') as f:
        X_test_seq = f['X_test_seq'][:]
        X_test_region = f['X_test_region'][:]
        y_test_RBP = f['y_test_RBP'][:]
        
        # 获取其他数据（如果存在）
        y_test_name = f['y_test_name'][:] if 'y_test_name' in f else None
        y_train = f['y_train'][:] if 'y_train' in f else None
    
    return X_test_seq, X_test_region, y_test_RBP, y_test_name, y_train

def get_top_predictions(predictions, labels, threshold=0.5, top_n=10):
    """
    获取高置信度的预测结果
    """
    positive_indices = np.where(labels == 1)[0]
    positive_predictions = predictions[positive_indices]
    
    # 筛选高于阈值的预测
    high_conf_mask = positive_predictions > threshold
    high_conf_indices = positive_indices[high_conf_mask]
    high_conf_scores = positive_predictions[high_conf_mask]
    
    # 按分数排序并获取top_n
    sorted_indices = np.argsort(-high_conf_scores)[:top_n]
    
    return high_conf_indices[sorted_indices], high_conf_scores[sorted_indices]
