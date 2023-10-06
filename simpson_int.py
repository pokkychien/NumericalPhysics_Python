"""
    simpson_int.py
    ~~~~~~~~~~~
    辛普森數值積分 (Simpson's rule) 的簡單範例。
    設一函數f(x)，上下界a與b，切割成n分，
    simpson(f,a,b,n) 將會對f(x)在a與b之間進行數值積分。
    
    此積分公式為：
    integral(f(x),a,b) = h/3 * (f(a)+f(b)+4*sum(f(a+2i*h))+2*sum(f(a+(2i+1)*h)))
    其中h=(b-a)/n，i=1,2,...,n-1。
    
    同時，exact_integral(a,b)為f(x)在a與b之間的解析解，將數值解與解析解相比較。
    例如說，f(x)=sin(x)在a與b之間積分的解析解為-cos(b)+cos(a)。

    Created by Chang Kai-Po @ Jian Lab, 2023/3/19
"""
import math

def f(x):
    """ 目標函數 """
    return math.sin(x)

def simpson (f,a,b,n):
    """ 辛普森積分的計算 """
    h = (b-a)/n
    s = f(a)+f(b) #邊界條件
    for i in range(1,n,2): #奇數項
        s += 4*f(a+i*h)
    for i in range(2,n-1,2): #偶數項
        s += 2*f(a+i*h)
    return s*h/3

def exact_integral(a,b):
    """ 解析解 """
    return -math.cos(b)+math.cos(a)

def main():
    a = 0.0
    b = math.pi
    n = 1000
    print("將 sin(x) 從 %g 到 %g 做積分" % (a,b))
    print("數值解: %g" % exact_integral(a,b))
    print("解析解: %g" % simpson(f,a,b,n))
    
if __name__ == "__main__":
    main()