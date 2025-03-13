import numpy as np

def matrice_iso_transversale(El, Et, Vlt, Vtt, Glt):
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

def matrice_orthotrope(E1, E2, E3, v12, v13, v23, G12, G23, G31):
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
    :param h: Épaisseur totale du stratifié
    :return: (A, B, D) matrices de rigidité
    """
    ply = len(angles)
   
    S = matrice_orthotrope_reduite(E1=E1, E2=E2, G12=G12, V12=V12)
     
    Q = np.linalg.inv(S)
    
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    
    j = 0
    k = 0
    
    for i in np.arange(h / ply, h + h / ply, h / ply):
        theta = np.deg2rad(angles[j])
        c = np.cos(theta)
        s = np.sin(theta)
        
        T = np.array([
            [c**2, s**2, -2*s*c],
            [s**2, c**2, 2*s*c],
            [s*c, -s*c, c**2 - s**2]
        ])
        
        T1 = np.array([
            [c**2, s**2, s*c],
            [s**2, c**2, -s*c],
            [-2*s*c, 2*s*c, c**2 - s**2]
        ])
        
        Q_global = T1 @ Q @ T
        
        A += (i - (i - h / ply)) * Q_global
        B += ((i**2) - (i - h / ply)**2) * Q_global
        D += ((i**3) - (i - h / ply)**3) * Q_global
        
        j += 1
        k += 1
    
    B *= 0.5
    D *= 1/3
    
    return A, B, D
