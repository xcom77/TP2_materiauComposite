import numpy as np


def auto_matrix(*args):
    """
    Sélectionne et appelle la fonction appropriée en fonction du nombre d'arguments donnés.

    :param args: Liste des paramètres donnés (peut être de taille variable)
    :return: (S, C) matrices de compliance et de rigidité
    """
    nb_args = len(args)

    if nb_args == 2:  # Matériau isotrope
        return matrice_isotrope(*args)
    
    elif nb_args == 5:  # Matériau isotrope transversal
        return matrice_iso_transversale(*args)
    
    elif nb_args == 9:  # Matériau orthotrope
        return matrice_orthotrope(*args)
    
    else:
        raise ValueError(f"Nombre d'arguments invalide : {nb_args}. "
                         "Attendu : 2 (isotrope), 5 (isotrope transversal), ou 9 (orthotrope).")


def matrice_iso_transversale(El, Et, Glt, Vlt, Vtt):
    """
    Calcule les matrices de compliance S et de rigidité C pour un matériau isotrope transversal.
    
    :param El: Module de Young longitudinal
    :param Et: Module de Young transverse
    :param Vlt: Coefficient de Poisson longitudinal-transverse
    :param Vtt: Coefficient de Poisson transverse-transverse
    :param Glt: Module de cisaillement longitudinal-transverse
    :return: (S, C) matrices de compliance et de rigidité
    """
    # Calcul de la matrice de compliance S
    S = np.array([
        [1/El, -Vlt/El, -Vlt/El, 0, 0, 0],
        [-Vlt/El, 1/Et, -Vtt/Et, 0, 0, 0],
        [-Vlt/El, -Vtt/Et, 1/Et, 0, 0, 0],
        [0, 0, 0, 2*(1+Vtt)/Et, 0, 0],
        [0, 0, 0, 0, 1/Glt, 0],
        [0, 0, 0, 0, 0, 1/Glt]
    ])
    
    # Calcul de la matrice de rigidité C
    C = np.linalg.inv(S)
    
    return S, C


def dec_trans_iso(matrice_C):
    """
    Calcule les propriétés mécaniques d'un matériau transiso en fonction de la matrice de rigidité C.
    
    :param matrice_C: numpy.ndarray, matrice de rigidité C
    :return: tuple (El, Et, Vtt, Vlt, Glt)
    """
    # Calcul de l'inverse de la matrice
    matrice = np.linalg.inv(matrice_C)
    
    # Tolérance pour comparer les valeurs flottantes
    tol = 1e-6 

    # Vérifications des conditions d'isotropie transversale
    if abs(matrice[4, 4] - matrice[5, 5]) > tol:
        raise ValueError("Erreur : les termes (5,5) et (6,6) sont incohérents.")
    if abs(matrice[1, 1] - matrice[2, 2]) > tol:
        raise ValueError("Erreur : les termes (2,2) et (3,3) sont incohérents.")
    if abs(matrice[2, 1] - matrice[1, 2]) > tol:
        raise ValueError("Erreur : les termes (3,2) et (2,3) sont incohérents.")
    if abs(matrice[0, 1] - matrice[0, 2]) > tol or abs(matrice[1, 0] - matrice[2, 0]) > tol or abs(matrice[0, 1] - matrice[1, 0]) > tol:
        raise ValueError("Erreur : les termes (1,2) et (1,3) sont incohérents.")

    # Calcul des propriétés mécaniques
    El = 1 / matrice[0, 0]
    Vlt = - matrice[1, 0] * El
    Et = 1 / matrice[1, 1]
    Vtt = - matrice[2, 1] * Et
    Glt = matrice[4, 4]

    # Vérification finale
    verif_value = 2 * ((1 + Vtt) / Et)
    if abs(matrice[3, 3] - verif_value) > tol:
        raise ValueError("Erreur : échec de la vérification des propriétés.")
    
    # Affichage des résultats
    print("Valeurs des coefficients de la matrice C :")
    print(f'La valeur de Et est : {Et:.5f}')
    print(f'La valeur de El est : {El:.5f}')
    print(f'La valeur de Vtt est : {Vtt:.5f}')
    print(f'La valeur de Vlt est : {Vlt:.5f}')
    print(f'La valeur de Glt est : {Glt:.5f}')
    
    return El, Et, Vtt, Vlt, Glt



def matrice_isotrope(E, v):
    """
    Calcule les matrices de compliance S et de rigidité C pour un matériau isotrope.
    
    :param E: Module de Young
    :param v: Coefficient de Poisson
    :return: (S, C) matrices de compliance et de rigidité
    """
    # Calcul de la matrice de compliance S
    S = np.array([
        [1/E, -v/E, -v/E, 0, 0, 0],
        [-v/E, 1/E, -v/E, 0, 0, 0],
        [-v/E, -v/E, 1/E, 0, 0, 0],
        [0, 0, 0, 2*(1+v)/E, 0, 0],
        [0, 0, 0, 0, 2*(1+v)/E, 0],
        [0, 0, 0, 0, 0, 2*(1+v)/E]
    ])
    
    # Calcul de la matrice de rigidité C
    C = np.linalg.inv(S)
    
    return S, C

def matrice_orthotrope(E1, E2, E3, G12, G23, G31, v12, v13, v23):
    """
    Calcule les matrices de compliance S et de rigidité C pour un matériau orthotrope.
    
    :param E1: Module de Young dans la direction 1
    :param E2: Module de Young dans la direction 2
    :param E3: Module de Young dans la direction 3
    :param v12: Coefficient de Poisson entre les directions 1 et 2
    :param v13: Coefficient de Poisson entre les directions 1 et 3
    :param v23: Coefficient de Poisson entre les directions 2 et 3
    :param G12: Module de cisaillement entre les directions 1 et 2
    :param G23: Module de cisaillement entre les directions 2 et 3
    :param G31: Module de cisaillement entre les directions 3 et 1
    :return: (S, C) matrices de compliance et de rigidité
    """
    # Compliance matrix S
    S = np.array([
        [1/E1, -v12/E1, -v13/E1, 0, 0, 0],
        [-v12/E1, 1/E2, -v23/E2, 0, 0, 0],
        [-v13/E1, -v23/E2, 1/E3, 0, 0, 0],
        [0, 0, 0, 1/G23, 0, 0],
        [0, 0, 0, 0, 1/G31, 0],
        [0, 0, 0, 0, 0, 1/G12]
    ])
    
    # Calcul de la matrice de rigidité C
    C = np.linalg.inv(S)
    
    return S, C


def matrice_iso_transversale_reduite(El, Et, Vlt, Glt):
    """
    Calcule les matrices de compliance S et de rigidité C réduites sous l'hypothèse de contrainte plane pour un matériau isotrope transversal.
    """
    S = np.array([
        [1/El, -Vlt/El, 0],
        [-Vlt/El, 1/Et, 0],
        [0, 0, 1/Glt]
    ])
    C = np.linalg.inv(S)
    return S, C

def matrice_isotrope_reduite(E, v):
    """
    Calcule les matrices de compliance S et de rigidité C réduites sous l'hypothèse de contrainte plane pour un matériau isotrope.
    """
    S = np.array([
        [1/E, -v/E, 0],
        [-v/E, 1/E, 0],
        [0, 0, 1/(2*(1+v)/E)]
    ])
    C = np.linalg.inv(S)
    return S, C

def matrice_orthotrope_reduite(E1, E2, v12, G12):
    """
    Calcule les matrices de compliance S et de rigidité C réduites sous l'hypothèse de contrainte plane pour un matériau orthotrope.
    """
    S = np.array([
        [1/E1, -v12/E1, 0],
        [-v12/E1, 1/E2, 0],
        [0, 0, 1/G12]
    ])
    C = np.linalg.inv(S)
    return S, C


def rotation_matrice_S(S_reduite, theta):
    """
    Effectue la transformation d'une matrice de compliance réduite S d'un repère local à un repère global.
    
    :param S_reduite: numpy.ndarray, Matrice de compliance réduite (3x3)
    :param theta: float, Angle de rotation en radians
    :return: numpy.ndarray, Matrice S dans le repère global
    """
    c = np.cos(theta)
    s = np.sin(theta)
    P = np.array([
        [c**2, s**2, -s*c],
        [s**2, c**2, s*c],
        [2*s*c, -2*s*c, c**2 - s**2]
    ])
    return P @ S_reduite @ P.T


def proprietes_mecaniques_apparentes(S_global):
    """
    Calcule les propriétés mécaniques apparentes à partir de la matrice de compliance globale.
    
    :param S_global: numpy.ndarray, Matrice de compliance globale (3x3)
    :return: tuple (E11, E22, nu12, G12)
    """
    E11 = 1 / S_global[0, 0]
    E22 = 1 / S_global[1, 1]
    nu12 = -S_global[0, 1] / S_global[0, 0]
    G12 = 1 / S_global[2, 2]
    return E11, E22, nu12, G12


def calcul_matrices_ABD(E1, E2, G12, V12, angles, h):
    """
    Calcule les matrices A, B et D pour un matériau stratifié en fonction des angles donnés.
    
    :param E1: Module de Young dans la direction 1 (GPa)
    :param E2: Module de Young dans la direction 2 (GPa)
    :param G12: Module de cisaillement (GPa)
    :param V12: Coefficient de Poisson
    :param angles: Liste des angles des plis en degrés
    :param h: Épaisseur d’un pli individuel
    :return: (A, B, D) matrices de rigidité arrondies à 4 décimales
    """
    ply = len(angles)  # Nombre de plis

    S, Q = matrice_orthotrope_reduite(E1=E1, E2=E2, G12=G12, v12=V12)
    Q = np.linalg.inv(S).astype(np.float64)

    A = np.zeros((3, 3), dtype=np.float64)
    B = np.zeros((3, 3), dtype=np.float64)
    D = np.zeros((3, 3), dtype=np.float64)

    for k in range(ply):
        theta = np.deg2rad(angles[k])  # Correction indexation
        
        # ⚠️ Correction de cos et sin pour 0° et 90° 
        if np.isclose(theta % (np.pi), 0):  # 0° ou 180°
            c, s = 1.0, 0.0
        elif np.isclose(theta % (np.pi), np.pi / 2):  # 90° ou 270°
            c, s = 0.0, 1.0
        else:
            c, s = np.cos(theta), np.sin(theta)

        T = np.array([
            [c**2, s**2, -2*s*c],
            [s**2, c**2, 2*s*c],
            [s*c, -s*c, c**2 - s**2]
        ], dtype=np.float64)
        
        T1 = np.array([
            [c**2, s**2, s*c],
            [s**2, c**2, -s*c],
            [-2*s*c, 2*s*c, c**2 - s**2]
        ], dtype=np.float64)
        
        # Produit matriciel avec np.dot() pour éviter les erreurs
        Q_global = np.dot(np.dot(T1, Q), T)

        Q_global = np.round(Q_global, 3)

        # Vérification des matrices pour les angles critiques
        #print(f"Angle: {np.rad2deg(theta)}°")
        #print(f"Q_global:\n{Q_global}\n")

        z_k = k * h  # Position inférieure du pli
        z_k1 = (k + 1) * h  # Position supérieure du pli

        A += (z_k1 - z_k) * Q_global
        B += (z_k1**2 - z_k**2) * Q_global
        D += (z_k1**3 - z_k**3) * Q_global

        # Arrondi après chaque itération
        A = np.round(A, 4)
        B = np.round(B, 4)
        D = np.round(D, 4)

    # Application des coefficients finaux APRÈS l’arrondi
    B *= 0.5
    D *= 1/3

    # Arrondi final
    A = np.round(A, 4)
    B = np.round(B, 4)
    D = np.round(D, 4)

    return A, B, D

def eshelby_tensor(nu_M):
    """
    Calcule la matrice du tenseur d'Eshelby pour une inclusion cylindrique
    dans un matériau isotrope caractérisé par le coefficient de Poisson nu_M.

    :param nu_M: Coefficient de Poisson du matériau isotrope
    :return: Matrice 6x6 du tenseur d'Eshelby
    """
    # Initialisation d'une matrice 6x6 remplie de zéros
    S = np.zeros((6, 6))

    # Calcul des valeurs à partir des formules données
    S[1, 1] = S[2, 2] = (5 - 4 * nu_M) / (8 * (1 - nu_M))
    S[1, 2] = S[2, 1] = (4 * nu_M - 1) / (8 * (1 - nu_M))
    S[1, 0] = S[2, 0] = nu_M / (2 * (1 - nu_M))
    S[3, 3] = (3 - 4 * nu_M) / (4 * (1 - nu_M))
    S[4, 4] = S[5, 5] = 1 / 2

    S = np.round(S,4)

    return S


def diluted_inclusion_localization(proprietes_mecanique_matrice,proprietes_mecanique_fibre):
    """
    Calcule le tenseur de localisation pour une inclusion diluée en utilisant la formule donnée.
    
    :param proprietes_mecanique_matrice: Paramètres mécaniques de la matrice (E, v, etc.)
    :param proprietes_mecanique_fibre: Paramètres mécaniques de la fibre (E, v, etc.)
    :return: Matrice A_dil du tenseur de localisation
    """
    # Calcul des matrices de compliance (S) et de rigidité (C) pour la matrice et la fibre
    S_m, C_m = auto_matrix(*proprietes_mecanique_matrice)  # Matrice
    S_l, C_l = auto_matrix(*proprietes_mecanique_fibre)    # Fibre

    # Calcul du tenseur d'Eshelby
    S_esh = eshelby_tensor(proprietes_mecanique_matrice[1])  

    # Calcul du terme (C^M)^-1
    S_m = np.linalg.inv(C_m)

    # Calcul du terme (C^I - C^M)
    delta_C = C_l - C_m

    # Calcul du tenseur de localisation A_dil
    A_dil = np.linalg.inv(np.eye(6) + np.einsum('ijkl,klmn,mnop->ijop', S_esh, S_m, delta_C))

    return A_dil



