"""
    local_extreme_bisect.py
    ~~~~~~~~~~~~~~~
    Finding local maxima and local minima of a function by bisection method, without using differential.    
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
    return x3, iter

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
    x, iter = local_extreme(f, a, b, tol, max_iter)
    print("在 x = %f, f(x) = %f 處找到極值，使用迴圈數為 %d。" % (x, f(x), iter))
          
    # Plot the function and the local maxima and local minima
    x_plot = np.linspace(a, b, 100)
    y_plot = f(x_plot)
    plt.title("Local maxima and local minima")
    plt.plot(x_plot, y_plot)
    plt.plot(x, f(x), 'ro')
    plt.show()

if __name__ == "__main__":
    main()