"""
    hw1_1014.py
    ~~~~~~~~~~~
    數值方法作業一，已知橢圓方程式為 x^2/4 + y^2/9 = 1，
    求函數f(x, y) = 2x-y的最大值。
    Created by Chang Kai-Po @ Jian Lab, NCTU, Taiwan, on 2023/3/26.
"""
import numpy as np
from matplotlib import pyplot as plt

def local_maxima_golden(f, a, b, tol, max_iter):
    """
    黃金比例法求最大值
    """
    # 設定初始點
    x_low, x_high = a, b
    phi = (1 + np.sqrt(5)) / 2 # 黃金分割比例
    x1 = x_low + (x_high - x_low) / phi # x_low和x_high之間的黃金分割點
    x2 = x_high + (x_low - x_high) / phi # x_high和x_low之間的黃金分割點
    f1, f2 = f(x1), f(x2) # 計算函數值
    iter = 0     
    li = [] 
    # 進行迭代
    while abs(x2-x1) > tol and iter < max_iter:          
        if f2 > f1: # f2的值較大
            x_high = x1
            x1 = x2 # x_low和x_high之間的黃金分割點
            x2 = x_high + (x_low - x_high) / phi # x_high和x_low之間的黃金分割點
        else: # f1的值較大
            x_low = x2
            x2 = x1 # x_high和x_low之間的黃金分割點
            x1 = x_low + (x_high - x_low) / phi # x_low和x_high之間的黃金分割點            
        f1, f2 = f(x1), f(x2) # 計算函數值
        iter += 1        
        li.append((x1, x2))
    return x_low, iter, li

def main():
    # 函數定義，推導後f(x, y) = 2x-y = 4*cos(t) - 3*sin(t)
    f = lambda t: 4*np.cos(t) - 3*np.sin(t)
    
    # 積分範圍與參數
    a, b, tol, max_iter = -5, 5, 1e-10, 100

    # 以黃金分割法找到函數的局部最大值
    t, iter, li = local_maxima_golden(f, a, b, tol, max_iter)
    print("在 t = %f, f(t) = %f 處找到極值，使用迴圈數為 %d。" % (t, f(t), iter))
    print("此時x=%f, y=%f" % (2*np.cos(t), 3*np.sin(t)))
    
    # 對x2的收斂情形及函數與極值進行繪圖
    x = [i[1] for i in li]    
    x_plot = np.linspace(a, b, 100)    
    y_plot = f(x_plot)    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].set_title("每次迭代的X2收斂情形")
    axs[1].set_title("函數與極值")
    axs[0].plot(x)       
    axs[1].plot(x_plot, y_plot)
    axs[1].plot(t, f(t), 'ro')
    fig.suptitle("黃金分割法找極值的效能")
    plt.rcParams['font.sans-serif'] = ['DFKai-SB'] # for Chinese characters
    plt.show()  
    
if __name__ == "__main__":
    main()