"""
    population.py
    ~~~~~~~~~~~~~
    dy/dt = ky - cy^2為具備外在競爭關係的的人口增長模型，
    本模組用來計算人口增長的數值解。

    This module contains a method for population growth estimation 
    modelled by an one-dimention differential equation, named logistic 
    equation of population growth.

    The differential equation is: dy/dt = ky - cy^2, and will be
    estimated by numerical methods. The initial condition is y(0) = 1.

    Chang Kai-Po @ Jian Lab 2023/03/13
"""
import math
import matplotlib.pyplot as plt

def finit_diff(y, k, c, h):
    """
    根據定義計算此微分方程的有限差分。
    若dy/dt = ky - cy^2，經過推導可得到
    f(x+h) = (1+hk)f(x) - hcf(x)^2
    """
    return (1+h*k)*y - h*c*y**2

def ranged_finit_diff (y0, xlow, xhigh, k, c, h):
    """
    透過有限差分計算此微分方程在在xlow與xhigh此一範圍之間的估計值。
    """
    y = [y0]
    level = int((xhigh - xlow)/h)
    for x in range(level):
        y.append(finit_diff(y[x], k, c, h))
    #print(y)        
    return y

def ranged_exact_sol(y0, xlow, xhigh, k, c, h):
    """
    傳回此微分方程精確解在xlow與xhigh此一範圍之間的結果。
    詳解請見 https://ch-hsieh.blogspot.com/2016/03/blog-post_10.html
    """
    y = [y0]
    level = int((xhigh - xlow)/h)
    for x in range(level):
        xvalue = xlow + x*h
        a = k*y0
        b = (k-(c*y0))*(math.exp((-k)*xvalue))
        d = c*y0
        y.append(a/(b+d))
    #print(y)        
    return y

def plot_finit_diff(y0, xlow, xhigh, k, c, h):
    """
    繪製有限差分的結果。
    """
    y = ranged_finit_diff(y0, xlow, xhigh, k, c, h)
    x = [xlow + i*h for i in range(len(y))]
    plt.plot(x, y)
    plt.show()

def plot_exact_sol(y0, xlow, xhigh, k, c, h):
    """
    繪製精確解的結果。
    """
    y= ranged_exact_sol(y0, xlow, xhigh, k, c, h)
    x = [xlow + i*h for i in range(len(y))]
    plt.plot(x, y)
    plt.show()

def plot_both_sol(y0, xlow, xhigh, k, c, h):
    """
    繪製有限差分與精確解的結果。
    """
    y1 = ranged_finit_diff(y0, xlow, xhigh, k, c, h)
    y2 = ranged_exact_sol(y0, xlow, xhigh, k, c, h)
    x = [xlow + i*h for i in range(len(y1))]
    plt.plot(x, y1, label='finit diff')
    plt.plot(x, y2, label='exact sol')
    plt.legend()
    plt.show()

#假設y(0)=10, c=0.1, h=0.1, xlow=0, xhigh=10
#plot_finit_diff(10, 0, 10, 0.1, 0.1, 0.0001) #繪製有限差分的結果
#plot_exact_sol(10, 0, 10, 0.1, 0.1, 0.0001) #繪製精確解的結果
plot_both_sol(1, 0, 50, 0.1, 0.01, 0.1) #繪製兩種解的結果