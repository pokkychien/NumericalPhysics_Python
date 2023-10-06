"""
heatconduction.py
~~~~~~~~~~~~~~~~~
加熱桿的熱傳導問題，使用有限差分法求解。
微分方程為：d^2T/dx^2 + k*(Ta-T) = 0
解析解為：T = 73.4523*exp(10kx) - 53.4523*exp(-10kx) + Ta
其中Ta為加熱桿本身的溫度，T為加熱桿的溫度分布。

Chang Kai-Po @ Jian Lab 2023/03/15
"""

import numpy as np
from scipy.sparse import csr_matrix 
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def coefficient_matrix (N, k, h, Ta):
    """
    建立線性方程式的係數矩陣，其中k為傳導係數，h為每次差分的差異，
    Ta為加熱棒。
    此差分方程式的解為 (2+kh^2) - f(x+h) - f(x-h) = kh^2Ta + c
    故係數矩陣為
    [2+kh^2, -1, 0, 0, ..., 0]
    [-1, 2+kh^2, -1, 0, ..., 0]
    [0, -1, 2+kh^2, -1, ..., 0]
    [...  ...  ...  ...  ...]
    [0, 0, ..., ..., -1, 2+kh^2]
    """
    A = np.zeros((N,N))
    A[0,0], A[0,1] = 2+k*h**2, -1
    A[N-1,N-1], A[N-1,N-2] = 2+k*h**2, -1
    for i in range(1,N-1):
        A[i,i-1] = -1
        A[i,i] = 2+k*h**2
        A[i,i+1] = -1
    return csr_matrix(A)

def heatconduction_numerical(L, T0, T1, Ta, k, N):
    """
    熱傳導問題的數值解，其中L為加熱桿的長度，T0為加熱桿左端的溫度，
    T1為加熱桿右端的溫度，Ta為加熱桿本身的溫度，k為傳導係數，
    N為加熱桿的分割數。
    """
    # 參數設定
    dx = L/(N+1) # 加熱桿的每一分割的長度    
    # 建立線性方程式的係數矩陣
    A = coefficient_matrix(N, k, dx, Ta)
    # 建立線性方程式的右邊向量
    b = np.repeat(k*dx**2*Ta, N)
    b[0] = T0
    b[N-1] = T1
    T = np.zeros(N+2)
    T[1:-1] = spsolve(A,b)
    T[0] = T0
    T[-1] = T1    
    return T    

def heatconduction_analytical(L, T0, T1, Ta, k, N):
    """
    熱傳導問題的數值解，其中L為加熱桿的長度，T0為加熱桿左端的溫度，
    T1為加熱桿右端的溫度，Ta為加熱桿本身的溫度，k為傳導係數，
    N為加熱桿的分割數。
    """
    x = np.linspace(0,L,N+2)
    y = 73.4523*np.exp(0.1*x) - 53.4523*np.exp(-0.1*x) + Ta #解析解            
    return y

if __name__ == '__main__':
    T0 = 40.0 # 左端的溫度
    T1 = 200.0 # 右端的溫度
    Ta = 20.0 # 加熱桿本身的溫度
    k = 0.01 # 熱傳導係數
    N = 4 # 加熱桿的分割數
    L = 10.0 # 加熱桿的長度

    # 比較數值解與解析解, N=4
    x1 = np.linspace(0,L,N+2)
    numerical = heatconduction_numerical(L,T0,T1,Ta,k,N)
    analytical = heatconduction_analytical(L,T0,T1,Ta,k,N)
    print("在N＝4時，數值解的{0}個x值分別為{1}".format(N, np.array_str(numerical[1:-1], precision=2, suppress_small=True)))
    print("在N＝4時，解析解的{0}個x值分別為{1}".format(N, np.array_str(analytical[1:-1], precision=2, suppress_small=True)))
    print("在N＝4時，數值解與解析解的誤差為{0}".format(np.array_str(numerical[1:-1]-analytical[1:-1], precision=2, suppress_small=True)))
    
    # 比較數值解與解析解, N=100
    N = 100
    x2 = np.linspace(0,L,N+2)
    numerical_2 = heatconduction_numerical(L,T0,T1,Ta,k,N)
    analytical_2 = heatconduction_analytical(L,T0,T1,Ta,k,N)
    print("在N＝100時，數值解的{0}個x值分別為{1}".format(N, np.array_str(numerical_2[1:-1], precision=4, suppress_small=True)))
    print("在N＝100時，解析解的{0}個x值分別為{1}".format(N, np.array_str(analytical_2[1:-1], precision=4, suppress_small=True)))
    print("在N＝100時，數值解與解析解的誤差為{0}".format(np.array_str(numerical_2[1:-1]-analytical_2[1:-1], precision=4, suppress_small=True)))
    
    # 繪圖      
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    # N=4時的數值解與解析解
    ax1.plot(x1, numerical, color='tab:blue', marker='o', markersize=2, label='Numerical')
    ax1.plot(x1, analytical, color='tab:red', marker='o', markersize=2, label='Analytical')
    ax1.set_title('N=4')

    # N=100時的數值解與解析解
    ax2.plot(x2, numerical_2, color='tab:blue', marker='o', markersize=2, label='Numerical')
    ax2.plot(x2, analytical_2, color='tab:red', marker='o', markersize=2, label='Analytical')
    ax2.set_title('N=100')
    
    # 設定圖例
    plt.xlabel('x')
    plt.ylabel('T')
    plt.show()    
