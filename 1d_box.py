"""
    1d_box.py
    ~~~~~~~~~
    Solve the particle in 1-dimentional box problem, and plot the wave function.
    Using numerical method for eigenvalue problem for Schrodinger equation.
"""
#modified by pokky
import numpy as np

def schrodinger (x):
    """
        Using numerical method for eigenvalue problem for Schrodinger equation.
    """
    # The solution matrix 
    M = np.zeros((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i == j:
                M[i][j] = 2
            elif i == j + 1 or i == j - 1:
                M[i][j] = -1