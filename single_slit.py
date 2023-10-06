"""
    single_slit.py
    ~~~~~~~~~~~~~~
    使用高斯積分進行單狹縫繞射計算。    
    Created by Chang Kai-Po @ Jian Lab, NCTU, Taiwan, on 2023/3/26.    
"""
import numpy as np 
from scipy.constants import milli, nano, c, pi
from scipy.integrate import quad #使用一般數值積分
from matplotlib import pyplot as plt

def slit_intensity (L, y, wavelength, d):
    """
    在泰勒展開近似後，單狹縫公式在某一個點的值為：
    I(y, x) = exp(i*2*pi*a*(y-x)^2)
    其中 a = k0 / 2L = 2*pi / (wavelength * 2L) = pi / (wavelength * L)  
    """
    #波長的單位為nm
    a = pi / (wavelength * nano * L) 
    #等等要對x做積分，所以先將y帶入公式
    f = lambda x: np.exp(1j*a*((y-x)**2))
    #積分範圍為[-d/2, d/2]    
    f_real = lambda x: np.real(f(x))
    f_imag = lambda x: np.imag(f(x))
    quad_real = quad(f_real, -d/2*milli, d/2*milli)[0]
    quad_imag = quad(f_imag, -d/2*milli, d/2*milli)[0]    
    return (abs(quad_real + 1j*quad_imag))**2

def main():
    y= np.linspace(-0.1, 0.1, 1000)
    z = np.zeros(1000)
    for i in range(1000):
        z[i] = slit_intensity(2, y[i], 630, 0.1)       
    plt.plot(y, z)
    plt.show()
    
if __name__ == '__main__':
    main()