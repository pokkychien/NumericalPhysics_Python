"""
    muller_root_finished.py
    Muller法找根(假設每個整數之間僅有一個根)
    註1: 如果兩個整數之間有超過一個根，那這支沒辦法解。
    註2: 這支不解虛根
    Chang Kai-Po @ Jian Lab 2023/03/08

"""
import math 
from scipy import constants

def test_case (x):  #用來測試的函數
    y = x**4 + x**3 - 2*x**2 + x - 6
    return y

def test_case_2 (x):  #用來測試的函數
    y = math.sin(3*x)
    return y

def muller_core (f, x0, x1, x2):
    """
    !!除了單元測試用途以外，請不要直接呼叫這個函數!!
    給予三個從小到大的x值: x0, x1, x2，並且f(x0), f(x1), f(x2)的符號不同
    算出下一個二次曲線的根
    """
    h1, h2 = (x1-x0), (x2-x1)
    delta1, delta2 = (f(x1)-f(x0))/h1, (f(x2)-f(x1))/h2
    d = (delta2 - delta1) / (h2 + h1)
    #以上是為了方便運算先設定的一些值
    if delta1 == delta2:                          #兩斜率相等，在同一直線上
        return x2 - ((x2-x1)/(f(x2)-f(x1))*f(x2)) #割線公式
    b = delta2 + h2*d
    D = (b**2 - 4*f(x2)*d)**0.5
    if D < 0:                                     #不討論虛根
        return;
    if abs(b-D) < abs(b+D):                       #選擇較大的值
        E = b + D
    else:
        E = b - D
    h = -2*f(x2)/E
    return x2 + h

def muller_root (f, x0, epsilon):
    """
    給予兩個從小到大的x值: x0, x2，並且f(x0)與f(x2)之間有根
    用Muller法找出根，精度為epsilon
    """
    x2 = x0 + 1
    x1 = x2 - ((x2-x0)/(f(x2)-f(x0))*f(x2))     #以割線做出第一個猜測值
    while f(x2)**2 > epsilon:                   #如果f(x2)的平方大於精度，則繼續
        x0, x1, x2 = x1, x2, muller_core (f, x0, x1, x2)
    return x2
    
def find_integer (f, low, high):
    """
      以low為下界，high為上界，搜尋函數f(x)在那些整數間可能會有根
      之後，傳回所有可能會有根的整數下緣
      (例如說如果有根在2-3之間與7-8之間則傳回[2,7])
    """
    li = []         #放一個空array填解
    for i in range(low, high):
        if f(i)*f(i+1) < 0:
            li.append(i);
    return li

def all_root (f, low, high, epsilon):
    """
      以low為下界, high為上界，搜尋函數f(x)在這個範圍內可能的根，
      精度為epsilon。      
    """
    introot = find_integer (f, low, high)      #搜尋函數f(x)在那些整數間可能會有根。
    #print(introot)
    li = []
    if not introot:
        return;                                #如果這範圍沒有，就結束
    for item in introot:
        li.append(muller_root (f, item, epsilon))
    return li

print(all_root(test_case, -20, 20, 0.0001))
print(all_root(test_case_2, -20, 20, 0.0001))
print(constants.G)