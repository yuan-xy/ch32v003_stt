import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体，确保中文正常显示
plt.rcParams["font.family"] = ["SimHei"]
plt.rcParams["axes.unicode_minus"] = False  # 解决负号显示问题

def hann_window(N):
    """生成汉宁窗函数"""
    n = np.arange(N)  # n的范围：0到N-1
    w = 0.5 * (1 - np.cos(2 * np.pi * n / (N - 1)))  # 汉宁窗公式
    return n, np.round(w*32).astype(np.int8) 

# 定义LPCNet窗口函数（用户提供的数组）
lpcnet_window = np.array([
    0,0,0,0,1,1,1,2,2,3,3,4,5,5,6,7,8,9,10,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,25,
    26,27,27,28,29,29,29,30,30,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
    32,32,32,32,32,32,32,32,32,32,32,32,31,31,31,31,31,30,30,29,29,29,28,27,27,26,25,25,24,23,22,
    21,20,19,18,17,16,15,14,13,12,11,10,10,9,8,7,6,5,5,4,3,3,2,2,1,1,1,0,0,0,0
])
N = 128  # 窗口长度，与数组长度一致
n = np.arange(N)  # 对应的n值（0-127）

# 生成汉宁窗数据
n_hann, w_hann = hann_window(N)
# breakpoint()

# 绘制对比图像
plt.figure(figsize=(12, 7))
# 绘制汉宁窗
plt.plot(n_hann, w_hann, 'b-', linewidth=2, label='汉宁窗 (Hann Window)')
# 绘制LPCNet窗口
plt.plot(n, lpcnet_window, 'r--', linewidth=2, label='LPCNet窗口函数')

# 图像标注
plt.title(f'汉宁窗与LPCNet窗口函数对比 (N=128)', fontsize=14)
plt.xlabel('n', fontsize=12)
plt.ylabel('归一化幅度', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0, N-1)
plt.ylim(0, 32.05)  # 最大值为1.0，留一点余量
plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
plt.legend(fontsize=12)  # 显示图例
plt.show()