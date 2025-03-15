from function import *

E1 = 132
E2 = 9.5
G12 = 5.2
V12 = 0.28

angle = [45,-45,-45,45]
h=0.256
A, B, D = calcul_matrices_ABD(E1=E1, E2=E2,G12=G12,V12=V12,angles=angle,h=h)

print("A1 : ",A)
print("B1 : ",B)
print("D1 : ",D)

angle = [45,45,45,45]
h=0.256
A, B, D = calcul_matrices_ABD(E1=E1, E2=E2,G12=G12,V12=V12,angles=angle,h=h)


print("A2 : ",A)
print("B2 : ",B)
print("D2 : ",D)

angle = [45,0,-45,90,90,-45,0,45]
h=0.256
A, B, D = calcul_matrices_ABD(E1=E1, E2=E2,G12=G12,V12=V12,angles=angle,h=h)


print("A3 : ",A)
print("B3 : ",B)
print("D3 : ",D)


nu_M_example = 0.3
S_matrix = eshelby_tensor(nu_M_example)
print("Tenseur d'Eshelby pour ν^M =", nu_M_example)
print(S_matrix)

C_m, C_l, Adil = diluted_inclusion_localization(proprietes_mecanique_fibre=[218,22,40,0.22,0.33],proprietes_mecanique_matrice=[3.6,0.39])

C_h = homogenized_rigidity_dilution(C_m, C_l, Adil,V_I=0.4)

print(C_h)
print(dec_trans_iso(C_h))


Vf = 0.4     # Fraction volumique de fibre (exemple)
Ef_L = 218    # GPa (Module de Young longitudinal des fibres)
Ef_T = 22     # GPa (Module de Young transverse des fibres)
Gf_LT = 40    # GPa (Module de cisaillement des fibres)
vf_LT = 0.22  # Coefficient de Poisson in-plane des fibres
vf_TT = 0.3   # Coefficient de Poisson transverse-transverse des fibres
Em = 3.6      # GPa (Module de Young de la matrice)
vm = 0.39     # Coefficient de Poisson de la matrice

proprietes_composite = rule_of_mixtures(Vf, Ef_L, Ef_T, Gf_LT, vf_LT, vf_TT, Em, vm)

for prop, valeur in proprietes_composite.items():
    print(f"{prop}: {valeur:.4f} GPa")
