"""
    simpson_vs_gauss.py
    ~~~~~~~~~~~~~~~~~~~
    比較辛普森積分和高斯積分的誤差。
    範例函數為 x^2*exp(x)。
    比較的範圍為 n= 3 到 10。
    
    Chang Kai-Po @ Jian Lab, NCTU, Taiwan, 2023/03/20
"""

import numpy as np
from matplotlib import pyplot as plt    

def gauss(f, a, b, n):
    """高斯積分"""
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * np.sum(w * f((b - a) / 2 * x + (b + a) / 2))

def simpson(f, a, b, n):
    """辛普森積分"""
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return h / 3 * (np.sum(y[0:-1:2]) + 4 * np.sum(y[1::2]) + np.sum(y[2::2]))

def exact(a, b):
    """解析解"""
    return (b**2 - 2*b + 2) * np.exp(b) - (a**2 - 2*a + 2) * np.exp(a)

def f(x):
    return x**2*np.exp(x)

def main():
    #計算在 n = 5 時的精確積分值，以及辛普森積分和高斯積分的值
    exact_integral = exact(0, 10)   
    simpson_integral = simpson(f, 0, 10, 10)
    gauss_integral = gauss(f, 0, 10, 5)
    
    #計算 n = 5 時的誤差
    simpson_error = abs(exact_integral - simpson_integral)/exact_integral
    gauss_error = abs(exact_integral - gauss_integral)/exact_integral
    
    #列印結果
    print("Value for exact integral of x^2*exp(x) from 0 to 10 at n = 5: %g" 
          % exact_integral)
    print("Value for exact integral of x^2*exp(x) from 0 to 10 at n = 5: %g" 
          % simpson_integral)
    print("Value for exact integral of x^2*exp(x) from 0 to 10 at n = 5: %g" 
          % gauss_integral)
    print("Simpson's rule error: %g" % simpson_error)
    print("Gaussian quadrature error: %g" % gauss_error)
    
    #計算n=3~10時的誤差
    simpson_errors = []
    gauss_errors = []
    for i in range(3, 10):
        simpson_errors.append(abs(exact_integral - simpson(f, 0, 10, 2*(i+1)))/exact_integral)
        gauss_errors.append(abs(exact_integral - gauss(f, 0, 10, (i+1)))/exact_integral)
    x = [i for i in range (3,10)]    
    #對誤差作圖    
    plt.plot(x, simpson_errors, label="Simpson's rule")    
    plt.plot(x, gauss_errors, label="Gaussian quadrature")
    plt.xlabel("n")
    plt.legend(loc='upper right')
    plt.show()

if __name__ == "__main__":
    main()