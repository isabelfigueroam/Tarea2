#!/usr/bin/env python3
import numpy as np

def gaussxw(N):
    """
    Calcula puntos de cuadratura y pesos de la regla de Gauss–Legendre en el intervalo estándar [-1, 1].

    Parameters
    ----------
    N : int
        Número de puntos (orden) de la cuadratura Gauss–Legendre.

    Returns
    -------
    x : numpy.ndarray
        Arreglo de forma (N,) con los puntos de colocación en [-1, 1].
    w : numpy.ndarray
        Arreglo de forma (N,) con los pesos asociados a cada punto.

    Examples
    --------
    >>> x, w = gaussxw(3)
    >>> x.shape, w.shape
    ((3,), (3,))
    """
    x, w = np.polynomial.legendre.leggauss(N)
    return x, w

def gaussxwab(a, b, x, w):
    """
    Escala puntos de cuadratura y pesos de Gauss–Legendre desde [-1, 1] a un intervalo [a, b].

    Parameters
    ----------
    a : float
        Extremo inferior del intervalo destino.
    b : float
        Extremo superior del intervalo destino.
    x : numpy.ndarray
        Puntos de colocación en el intervalo estándar [-1, 1].
    w : numpy.ndarray
        Pesos correspondientes en el intervalo estándar [-1, 1].

    Returns
    -------
    x_scaled : numpy.ndarray
        Puntos de colocación escalados al intervalo [a, b].
    w_scaled : numpy.ndarray
        Pesos ajustados para el intervalo [a, b].

    Examples
    --------
    >>> x, w = gaussxw(2)
    >>> xs, ws = gaussxwab(1.0, 3.0, x, w)
    >>> (xs.min() >= 1.0) and (xs.max() <= 3.0)
    True
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def funcIntegrando(x):
    r"""
    Función integrando de la tarea:

    .. math::

        f(x) = x^6 - x^2 \sin(2x).

    Parameters
    ----------
    x : float or numpy.ndarray
        Punto(s) donde evaluar el integrando.

    Returns
    -------
    float or numpy.ndarray
        Valor(es) de :math:`f(x)`.

    Examples
    --------
    >>> out = funcIntegrando(1.0)
    >>> isinstance(out, float)
    True
    >>> arr = funcIntegrando(np.array([1.0, 2.0]))
    >>> arr.shape
    (2,)
    """
    return x**6 - (x**2 * np.sin(2 * x))

def integrar_gauss(a, b, N):
    """
    Integra `funcIntegrando(x)` en el intervalo [a, b] usando cuadratura Gauss–Legendre.

    Parameters
    ----------
    a : float
        Límite inferior de integración.
    b : float
        Límite superior de integración.
    N : int
        Número de puntos (orden) de la cuadratura.

    Returns
    -------
    float
        Aproximación de la integral de `funcIntegrando` en [a, b] con N puntos.

    Examples
    --------
    >>> val4 = integrar_gauss(1.0, 3.0, 4)
    >>> val6 = integrar_gauss(1.0, 3.0, 6)
    >>> abs(val6 - val4) < 1e-2
    True
    """
    x, w = gaussxw(N)
    xs, ws = gaussxwab(a, b, x, w)
    return float(np.sum(funcIntegrando(xs) * ws))
