"""
    least_square.py
    ~~~~~~~~~~~~~~
    用最小平方法對特定資料點進行二次函數近似，並繪製出來。
    Using least square method to approximate the function as a parabolic function.
    Created by Chang Kai-Po @ Jian Lab, NCTU, Taiwan, on 2023/4/6.    
"""
import numpy as np 
from matplotlib import pyplot as plt

def main():
    # 資料點，來自https://zh.wikipedia.org/zh-tw/%E6%9C%80%E5%B0%8F%E4%BA%8C%E4%B9%98%E6%B3%95
    x = np.array([208, 152, 113, 227, 137, 238, 178, 104, 191, 130])
    y = np.array([21.6, 15.5, 10.4, 31.0, 13.0, 32.4, 19.0, 10.4, 19.0, 11.8])
    
    # 排序
    x = np.sort(x)
    y = np.sort(y)
    
    # 計算最小平方法的係數
    A = np.vstack([x**2, x, np.ones(len(x))])
    # 通解為a = (A^T*A)^-1*A^T*y
    # 所以先計算X^T*X的逆矩陣
    B = np.linalg.inv(np.dot(A, A.T))
    # 最後計算通解
    a = np.dot(B, np.dot(A, y))
    m, c, b = a
    print("二次曲線近似為y=", m, "* x^2 + ", c, "* x + ", b, sep="")

    # 繪製圖形
    plt.plot(x, y, 'o', label='Original data', markersize=10)
    plt.plot(x, m*x**2 + c*x + b, 'r', label='Fitted line')
    plt.legend()
    plt.show()
    
if __name__ == '__main__':
    main()