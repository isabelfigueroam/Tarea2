#!/usr/bin/env python3

from scipy.special import legendre
import numpy as np

def gaussxw(N):
    """
    Calcula los puntos y pesos de la cuadratura de Gauss-Legendre
    en el intervalo estándar [-1, 1].

    Args:
        N (int): Número de puntos de cuadratura.

    Returns:
        tuple: (x, w)
            - x (numpy.ndarray): Puntos de colocación.
            - w (numpy.ndarray): Pesos asociados.

    Example:
        >>> x, w = gaussxw(3)
        >>> len(x), len(w)
        (3, 3)
    """
    x, w = np.polynomial.legendre.leggauss(N)
    return x, w


def gaussxwab(a, b, x, w):
    """
    Escala los puntos y pesos de la cuadratura de Gauss-Legendre
    desde [-1, 1] hasta un intervalo arbitrario [a, b].

    Args:
        a (float): Extremo inferior del intervalo.
        b (float): Extremo superior del intervalo.
        x (numpy.ndarray): Puntos en [-1, 1].
        w (numpy.ndarray): Pesos en [-1, 1].

    Returns:
        tuple: (x_scaled, w_scaled)
            - x_scaled (numpy.ndarray): Puntos escalados a [a, b].
            - w_scaled (numpy.ndarray): Pesos ajustados a [a, b].

    Example:
        >>> x, w = gaussxw(2)
        >>> xs, ws = gaussxwab(1, 3, x, w)
        >>> len(xs), len(ws)
        (2, 2)
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def funcIntegrando(x):
    """
    Define la función integrando para el problema planteado:

        f(x) = x^6 - x^2 * sin(2x)

    Args:
        x (float or numpy.ndarray): Valor o arreglo de valores donde se evalúa f(x).

    Returns:
        float or numpy.ndarray: Resultado de evaluar la función.

    Example:
        >>> funcIntegrando(1.0)
        0.09070257317431829
    """
    return x**6 - (x**2 * np.sin(2 * x))


def integrar_gauss(a, b, N):
    """
    Calcula la integral de la función definida en `funcIntegrando(x)`
    usando cuadratura de Gauss-Legendre con N puntos.

    Args:
        a (float): Límite inferior de integración.
        b (float): Límite superior de integración.
        N (int): Número de puntos de cuadratura.

    Returns:
        float: Aproximación de la integral en [a, b].

    Example:
        >>> integrar_gauss(1, 3, 4)
        2.639399999999999  # resultado aproximado

    Notes:
        - Para probar la convergencia, puede compararse este valor con el obtenido
          usando un N mucho mayor (ej. N=20).  
        - A medida que se incrementa N, la aproximación converge al valor de referencia.  
        - El valor “correcto” en este contexto se identifica como aquel en que
          el resultado deja de cambiar en las cifras significativas.
    """
    x, w = gaussxw(N)
    xs, ws = gaussxwab(a, b, x, w)
    return np.sum(funcIntegrando(xs) * ws)

