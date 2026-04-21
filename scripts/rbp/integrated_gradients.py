"""
简化的集成梯度实现
"""
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model

class IntegratedGradients:
    def __init__(self, model, steps=50):
        self.model = model
        self.steps = steps
    
    def explain(self, inputs, outc=0, reference=True, baseline=None):
        """
        计算集成梯度
        inputs: 输入数据
        outc: 输出类别索引
        reference: 是否使用零基准
        baseline: 自定义基准
        """
        # 获取模型输入
        if isinstance(inputs, list):
            input_tensors = [tf.convert_to_tensor(inp, dtype=tf.float32) 
                           for inp in inputs]
        else:
            input_tensors = [tf.convert_to_tensor(inputs, dtype=tf.float32)]
        
        # 创建基准
        if baseline is None:
            if reference:
                baseline = [tf.zeros_like(tensor) for tensor in input_tensors]
            else:
                baseline = [0.25 * tf.ones_like(tensor) for tensor in input_tensors]
        
        # 计算梯度
        attributions = self._compute_attributions(input_tensors, baseline, outc)
        
        return attributions
    
    def _compute_attributions(self, inputs, baseline, outc):
        """计算归因值"""
        # 简化的实现 - 可以根据需要扩展
        gradients = []
        
        for alpha in np.linspace(0, 1, self.steps):
            # 插值输入
            interpolated = []
            for inp, base in zip(inputs, baseline):
                interpolated.append(base + alpha * (inp - base))
            
            # 计算梯度
            with tf.GradientTape() as tape:
                for inp in interpolated:
                    tape.watch(inp)
                predictions = self.model(interpolated)
                output = predictions[:, outc]
            
            grad = tape.gradient(output, interpolated)
            if grad[0] is not None:
                gradients.append(grad[0].numpy())
        
        # 平均梯度
        if gradients:
            avg_grad = np.mean(gradients, axis=0)
            attribution = (inputs[0].numpy() - baseline[0].numpy()) * avg_grad
        else:
            attribution = np.zeros_like(inputs[0].numpy())
        
        return attribution
