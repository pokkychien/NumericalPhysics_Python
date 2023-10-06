"""
    potential_well.py
    ~~~~~~~~~~~~~~~~
    使用有限差分法估計一個電子在一階位能井中的波函數機率分布，
    以及其對應的能量。
    位能井的長度為 1，電子的質量為electron_mass，位能井的深度為無限。

    Created by Chang Kai-Po @ Jian Lab, 2023/03/18
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
from scipy.constants import electron_mass, hbar
import pyarma as pa 

def coeff_matrix(N):
    """
    有限差分法的係數矩陣。
    方程式為-(hbar^2 / 2m) d^2(psi(x) / dx^2) = E psi(x), 0 < x < L
    邊界條件為psi(0) = psi(L) = 0
    運算後，有限差分法的係數矩陣為
    A = [2 -1 0 0 0 ... 0]
        [-1 2 -1 0 0 ... 0]
        [0 -1 2 -1 0 ... 0]
        [.................]       
        [0 0 ... -1 2 -1]
        [0 0 ... ... -1 2]
    因此為一特徵值問題，將其轉換為求解A的特徵值和特徵向量。
    其中特徵值為E*2*epsilon^2*m/hbar^2，特徵向量為psi(x)。
    """
    A = np.zeros((N, N))
    A[0, 0] = 2
    A[N-1, N-1] = 2
    A[0, 1] = -1
    A[N-1, N-2] = -1
    for i in range(1, N-1):
        A[i, i-1] = -1
        A[i, i] = 2
        A[i, i+1] = -1 
    eigenvalues, eigenvectors =  pa.eig_gen(A, k=N-2, which='SM')
    sorted_indices = np.argsort(eigenvalues)    
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    return sorted_eigenvalues, sorted_eigenvectors     


eigenvalues = coeff_matrix(500)
epsilon = 0.002
print("First ten eigenvalues:", [f"{np.real(e):.7f}" for e in eigenvalues[0][:10]])

fig, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
for i in range(3):
    # Energy = eigenvalue *0.5* hbar^2 / epsilon^2 / m
    energy = eigenvalues[0][i] * 0.5 * hbar**2 / epsilon**2 / electron_mass
    axs[i].plot(eigenvalues[1][i], label=f"n={i+1} E={np.real(energy)}")
    axs[i].set_ylabel(r"$\psi(x)$")
    axs[i].legend(loc='upper right')
axs[2].set_xlabel(r"$x$")
plt.show()