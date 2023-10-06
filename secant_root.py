"""
    secant_root.py
    割線法找根(假設每個整數之間僅有一個根)
    註: 如果兩個整數之間有超過一個根，那這支沒辦法解。
    Chang Kai-Po @ Jian Lab 2023/03/03
"""

def test_case (x):  #用來測試的函數
    y = x**4 + x**3 - 2*x**2 + x - 6
    return y

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

def one_root (f, i, epsilon):
    """
      !!除了單元測試用途以外，請不要直接呼叫這個函數!!
      在已經知道函數f在i與i+1之間有一個根的前提下，
      用二分法逼近出這個根在哪裡，精度為epsilon。    
    """
    low, mid, high = float(i), float(i)+0.5, float(i)+1
    while f(mid)**2 > epsilon:
        mid = high - ((high-low)/(f(high)-f(low))*f(high)) #割線公式
        if f(low)*f(mid) < 0:          #根在下半，故異號
            high = mid                 #中間值變成上界            
        else:                          #根在上半，故異號
            low = mid;                 #中間值變成下界
    return mid                         #以下界為根 

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
        li.append(one_root (f, item, epsilon))
    return li

print("函數 x^4 + x^3 - 2*x^2 + x - 6 \n 在-20與20之間的根可能有:\n");
print(all_root(test_case, -20, 20, 0.0001))