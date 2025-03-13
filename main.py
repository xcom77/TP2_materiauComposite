from function import *

E1 = 132
E2 = 9.5
G12 = 5.2
V12 = 0.28
angle = [45,-45,-45,45]
h=0.256
A, B, D = calcul_matrices_ABD(E1=E1, E2=E2,G12=G12,V12=V12,angles=angle,h=h)

print(A)
print(B)
print(D)