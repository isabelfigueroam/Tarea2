# Tutorial de uso

Este módulo permite calcular la integral

\[
\int_{1}^{3} \left(x^6 - x^2\sin(2x)\right)\,dx
\]

usando cuadratura Gauss-Legendre. A continuación se muestra un ejemplo completo de uso.

---

## Ejemplo 

Este ejemplo recorre valores de N hasta que dos aproximaciones consecutivas difieren menos de `1e-6` (el valor de tolerancia).

```python
from cuadrature import integrar_gauss

a, b = 1.0, 3.0
tol = 1e-6

N = 1
prev = None

while True:
I_N = integrar_gauss(a, b, N)
if prev is not None:
diff = abs(I_N - prev)
print(f"N={N:2d} | I={I_N:.12f} | |Δ|={diff:.2e}")
if diff < tol:
print(f"\nConvergencia alcanzada con N = {N}")
break
else:
print(f"N={N:2d} | I={I_N:.12f} | |Δ|=n/a")
prev = I_N
N += 1
```

**Salida:**
```
N= 1 | I=134.054419962463 | |Δ|=n/a
N= 2 | I=306.819934495920 | |Δ|=1.73e+02
N= 3 | I=317.264151733829 | |Δ|=1.04e+01
N= 4 | I=317.345390334158 | |Δ|=8.24e-03
N= 5 | I=317.344226721970 | |Δ|=1.16e-03
N= 6 | I=317.344246889996 | |Δ|=2.02e-05
N= 7 | I=317.344246672227 | |Δ|=2.18e-07

Convergencia alcanzada con N = 7
```

---


## ¿Qué tan grande debe ser N?

En Gauss–Legendre, un orden N integra exactamente polinomios de grado ≤ (2N − 1).
Aquí el integrando tiene una parte polinómica `$x^6$` y una parte no polinómica `$−x^2·sin(2x)$`. Por eso el error no se anula con un N finito, pero decrece rápidamente al aumentar N. Entonces, buscamos un N tal que cumpla la tolerancia requerida por el problema.

---
