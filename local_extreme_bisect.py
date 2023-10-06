"""
    local_extreme_bisect.py
    ~~~~~~~~~~~~~~~
    以二分法找出函數的局部最大值和局部最小值。
    範例函數為 f(x) = -(x-3)**2+1。
    
    Created by Chang Kai-Po @ Jian Lab, NCTU, Taiwan, on 2023/3/21.
"""
import numpy as np
import matplotlib.pyplot as plt

def local_extreme(f, a, b, tol, max_iter):
    """
    Find local maxima and local minima of a function f by bisection method.
    """
    # 設定初始點
    x1, x2 = a, b
    x3 = (x1 + x2) / 2
    f1, f2, f3 = f(x1), f(x2), f(x3)
    iter = 0  
    li = []  
    # 進行迭代
    while abs(x2 - x1) > tol and iter < max_iter:    
        x4 = (x1 + x3) / 2 # x1和x3之間的中心點        
        x5 = (x3 + x2) / 2 # x2和x3之間的中心點
        f4, f5 = f(x4), f(x5)
        left_exam = (f4-f1)*(f4-f3) #左邊極值檢核
        right_exam = (f5-f3)*(f5-f2) #右邊極值檢核
        if left_exam * right_exam > 0: #同正或同負時，選擇值較大的一邊
            if left_exam > right_exam: #左邊可能性大於右邊
                x2, f2 = x3, f3
            else:   
                x1, f1 = x3, f3
        elif left_exam > 0: #只有可能在左邊
            x2, f2 = x3, f3 
        else: #只有可能在右邊
            x1, f1 = x3, f3
        x3 = (x1 + x2) / 2        
        f3 = f(x3)       
        iter += 1
        li.append(x3)
    return x3, iter, li

def main():
    # Define the function
    f = lambda x: -(x-3)**2+1

    # Define the interval
    a = -20
    b = 20

    # Define the tolerance and the maximum number of iterations
    tol = 1e-10
    max_iter = 100

    # Find the local maxima and local minima
    x, iter, li = local_extreme(f, a, b, tol, max_iter)
    print("在 x = %f, f(x) = %f 處找到極值，使用迴圈數為 %d。" % (x, f(x), iter))
    
    #對x3的收斂情形及函數與極值進行繪圖      
    x_plot = np.linspace(a, b, 100)
    y_plot = f(x_plot)
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].set_title("每次迭代的X3收斂情形")
    axs[1].set_title("函數與極值")
    axs[0].plot(li)
    axs[1].plot(x_plot, y_plot)
    axs[1].plot(x, f(x), 'ro')
    fig.suptitle("二分法找極值的效能")
    plt.rcParams['font.sans-serif'] = ['DFKai-SB'] # for Chinese characters
    plt.show()  

if __name__ == "__main__":
    main()