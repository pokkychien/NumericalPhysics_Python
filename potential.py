"""
potential.py
~~~~~~~~~~~~~~~~~
Finite differential solution for one-dimention infinite potential well.
Author: Chang Kai-Po @ Jian Lab 2023/03/16

"""

import numpy as np
from scipy.sparse import csr_matrix 
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def coefficient_well (N, h, V):
    """
    建立線性方程式的係數矩陣，其中V為無限深势能井的势能，
    h為每次差分的差異。
    此差分方程式的解為 (2+Vh^2) - f(x+h) - f(x-h) = 0
    故係數矩陣為
    [2+Vh^2, -1, 0, 0, ..., 0]
    [-1, 2+Vh^2, -1, 0, ..., 0]
    [0, -1, 2+Vh^2, -1, ..., 0]
    [...  ...  ...  ...  ...]
    [0, 0, ..., ..., -1, 2+Vh^2]
    """
    A = np.zeros((N,N))
    A[0,0], A[0,1] = 2+V*h**2, -1
    A[N-1,N-1], A[N-1,N-2] = 2+V*h**2, -1
    for i in range(1,N-1):
        A[i,i-1] = -1
        A[i,i] = 2+V*h**2
        A[i,i+1] = -1
    return csr_matrix(A)