# simple-kalman

from simple_kalman import Kalman

A = [[0,       1],
     [-k/m, -b/m]]
B = [[0],
     [1]]
C = [[1, 0]]    
D = 0

Kalman(A, B, C, D)
