# Método Numérico: Cuadratura Gauss-Legendre

La cuadratura Gaussiana es un método numérico para aproximar integrales definidas de la forma:

\[
\int_{a}^{b} f(x)\,dx \approx \sum_{i=1}^{N} w_i\, f(x_i)
\]

donde:

- \(x_i\) son los puntos de colocación.
- \(w_i\) son los pesos asociados a cada punto.
- \(N\) es el orden de la cuadratura (la cantidad de puntos utilizados).

A diferencia de otros métodos vistos en el curso como trapezoides o Simpson, que usan puntos equiespaciados, la cuadratura gaussiana elige los puntos \(x_i\) de manera mas óptima, para maximizar la precisión con el menor número de evaluaciones de la función.

En el caso de Gauss-Legendre, los puntos \(x_i\) se eligen como las raíces del polinomio de Legendre \(P_N(x)\), definido por la relación de recurrencia:

\[
(n+1)P_{n+1}(x) = (2n+1)x\,P_n(x) - n\,P_{n-1}(x), \quad P_0(x)=1,\quad P_1(x)=x
\]

Los pesos correspondientes están dados por:

\[
w_i = \frac{2}{\left(1 - x_i^2\right)\left[P_N'(x_i)\right]^2}
\]

La cuadratura Gauss-Legendre de orden \(N\) integra exactamente cualquier polinomio de grado ≤ \(2N - 1\).

Por ejemplo, con N = 3, se puede integrar exactamente cualquier polinomio de hasta grado 5, sin error.

Este método tiene la ventaja de ser rápido, preciso y sistemático. Con solo aumentar \(N\) mejora la aproximación hasta que el resultado se estabiliza dentro de una tolerancia requerida. En la práctica, se calcula con varios valores de N y se elige el primero que alcanza la convergencia (como se muestra en la sección *Tutorial de uso*).
