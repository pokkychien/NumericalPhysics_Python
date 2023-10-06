"""
    potential_well_out_of_box.py
    ~~~~~~~~~~~~~~~~
    使用有限差分法估計一個電子在一階位能井中的波函數機率分布，
    以及其對應的能量。
    位能井的長度為 1，電子的質量為electron_mass，位能井的深度為V0。

    Created by Chang Kai-Po @ Jian Lab, 2023/03/18
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import electron_mass, hbar

def coeff_matrix_out(N, boundary=10, V0=0.5):
    """
    有限差分法的係數矩陣。
    方程式為
    方塊內: -(hbar^2 / 2m) d^2(psi(x) / dx^2) = E psi(x), -L/2 < x < L/2
    方塊外: -(hbar^2 / 2m) d^2(psi(x) / dx^2) = E psi(x), -boundary*L < x < +boundary*L
    邊界條件為psi(-boundary*L) = psi(+boundary*L) = 0
    運算後，有限差分法的係數矩陣為
    A = [2 -1 0 0 0 ... 0]
        [-1 2 -1 0 0 ... 0]
        [0 -1 2 -1 0 ... 0]
        [.................]       
        [0 0 ... -1 2 -1]
        [0 0 ... ... -1 2]
    當x在井內時V(x) = 0，當x在井外時V(x) = V0
    因此為一特徵值問題，將其轉換為求解A的特徵值和特徵向量。
    其中特徵值為E*2*epsilon^2*m/hbar^2，特徵向量為psi(x)。
    """    
    A = np.zeros((N, N))
    A[0, 0] = 2+V0
    A[N-1, N-1] = 2+V0
    A[0, 1] = -1
    A[N-1, N-2] = -1
    for i in range(1, N-1):
        A[i, i-1] = -1
        if abs(i-N/2) < (N/(4*boundary+2)): #在井內
            A[i, i] = 2
        else:                      #在井外
            A[i, i] = 2+V0
        A[i, i+1] = -1    
    eigenvalues, eigenvectors = np.linalg.eig(A)
    sorted_indices = np.argsort(eigenvalues)    
    sorted_eigenvalues = eigenvalues[sorted_indices]    
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    # Get the absolute values of the eigenvectors
    sorted_eigenvectors = np.abs(sorted_eigenvectors)
    # Get the absolute values of the eigenvalues
    sorted_eigenvalues = np.abs(sorted_eigenvalues)
    return sorted_eigenvalues, sorted_eigenvectors   

eigenvalues = coeff_matrix_out(200, boundary=10, V0=1)
epsilon = 0.002
print("First ten eigenvalues:", [f"{np.real(e):.7f}" for e in eigenvalues[0][:10]])
print(eigenvalues[1][:10])

fig, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
for i in range(3):
    # Energy = eigenvalue *0.5* hbar^2 / epsilon^2 / m
    energy = eigenvalues[0][i] * 0.5 * hbar**2 / epsilon**2 / electron_mass
    axs[i].plot(eigenvalues[1][i], label=f"n={i+1} E={np.real(energy)}")
    axs[i].set_ylabel(r"$\psi(x)$")
    axs[i].legend(loc='upper right')    
axs[2].set_xlabel(r"$x$")
plt.show()