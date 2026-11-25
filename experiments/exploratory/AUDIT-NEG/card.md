# AUDIT-NEG: auditoría retroactiva de los negativos registrados (regla 39)

Ficha escrita el 24 de noviembre de 2026, **antes** de correr ningún control. Pedido del dueño (2026-11-23: "prueba
todo antes de darlo como negativo y verifica errores de implementacion"; "SIII HAZLOOOO" a la propuesta de auditar
todos los negativos). Hilo "Auditar los negativos por errores de implementación". Exploratorio: no cambia ningún
número registrado ni relaja ningún criterio (reglas 14, 15, 36). Un negativo registrado sigue registrado; si un
control encuentra un defecto, se anota "negativo con un defecto de implementación encontrado después" con fecha.

## Alcance

| línea | ficha | código principal |
|---|---|---|
| M5-ACAP (G-M5a, G-M5b) | `M5-ACAP/card.md` | `run_a.py` (Task 2M), `run_b.py` (oráculo C18c sobre `r0_run.Net`) |
| M5-EVALUE (G-M5c, 3 intentos) | `M5-EVALUE/card.md` | `src/ttmt/v2/m5.py`, `run_c.py` |
| M2-GP (G-M2) | `M2-GP/card.md` | `src/ttmt/v2/m2.py`, `run.py` |
| N-1 (N1-STACK) | `N1-STACK/card.md` | `run.py` |
| N-2 (N2-REPICK) | `N2-REPICK/card.md` | `run.py` |
| R0, R1, C18a | `results/taskR0`, `results/taskR1`, E50 | `code/r0`, `code/r1`, `src/ttmt/r1.py` |
| M1A-ONSET | `M1A-ONSET/card.md` | no es un negativo (pasó en EVAL); se revisa solo el blanqueo |

Fuera de alcance: M7 (aplazado por el dueño) y la línea M6-SPLIT (la audita su hilo).

Patrones de defecto que se buscan, todos ya confirmados en otra línea:
1. **Blanqueo** (PR #68, ADR-0012): datos divididos por la sd de cada estación y subespacio Q construido con G sin
   pesos. Se busca `station_scale` / división por sd por estación junto a `radiation_operator`/`basis` sin pesos.
2. **Retardo por interpolación lineal** (generador de M6-SPLIT): `np.interp` para un retardo fraccionario, que
   pierde energía en la componente retrasada.
3. **Búsqueda local que arranca de un solo punto** (`free_one_event` de M6-SPLIT).
4. **Operador de rayo recto** (FE-FIX): medio homogéneo en lugar de un modelo 1-D.

## Controles en la nube (sin portátil, sin datos del proyecto, datos sintéticos con semillas fijas)

Guion: `controls.py` (una función por control; salida `controls_out.txt`). Semilla base 20261124. Si un control
positivo o negativo no sale como se espera, la lectura de esa línea es "sin auditar: el control falla" y no se
concluye nada más sobre ella (regla 40).

### C-M5. El e-valor con la verdad conocida (`ttmt.v2.m5.evalue_mix`)
Las pruebas de `tests/test_m5.py` comprueban el nivel solo a p = 1e-1 y 1e-2; el gate G-M5c es a 1e-3.
K = 16 bloques (8 estaciones x P, S), operadores A ~ N(0, 0.3^2), J = 500 tensores de `prior_draws` al nivel del
ruido (como el test existente), eps = 0.05, s0 = 1, gamma = 0.3, H = 200 000 hipótesis.
- **M5-a (positivo de validez):** lecturas N(0, I_3) y fondo gaussiano exacto. Tasa de S >= log(1/p) para p = 1e-2 y
  1e-3, producto y promedio por estación (grupos de 2 bloques). Criterio: las cuatro tasas <= p, y la media de
  exp(S_avg) en [0.9, 1.1].
- **M5-b (negativo):** las mismas lecturas con un fondo demasiado ligero (gaussiano de sd 0.8). Criterio: la tasa del
  producto a 1e-3 > 1e-3 (la prueba sabe detectar un fondo malo).
- **M5-c (positivo de potencia):** lecturas = A m + N(0, I) con m del mismo pi multiplicado por 3. Criterio: la
  fracción con S_avg >= log(1000) es >= 0.5.
- **M5-d (descriptivo; decide la lectura del tercer intento):** lecturas de pico (máximo de energía de 9 muestras
  gaussianas, el lector registrado en forma ideal), fondo `RadialBackground.fit(left="flat")` ajustado sobre todas
  las lecturas (`--weighted`) de una muestra independiente de tamaño n_fit = 30 000 (orden del ajuste real por
  estación) y n_fit = 2 000 000 (casi exacto). Tasa de S_avg y del producto a 1e-3.
  Lectura fijada: si con n_fit grande la tasa es <= 1e-3 y con n_fit 30 000 es > 1e-3, el fallo de S_avg viene del
  error finito del fondo estimado (el negativo se sostiene como propiedad del método con fondo estimado); si con
  n_fit grande también es > 1e-3, hay un defecto en `RadialBackground` o en `evalue_mix` (contaminado) y se busca.

### C-M2. El GP con un campo conocido (`ttmt.v2.m2.GP1`, `Const1`, `time_blocks`)
Geometría sintética parecida a Ridgecrest: 3 000 eventos en una falla de 60 x 3 km, 2-12 km de profundidad, en orden
de tiempo a lo largo de la falla con mezcla (tiempos al azar); 6 pares (estación, fase). Residuo = constante +
campo Matern 3/2 (ell 8 km, sf 0.15 s) + ruido blanco (sn 0.15 s). Mismos 5 bloques de tiempo y mismo
`summ` que `M2-GP/run.py` (MSE drop, cobertura +-2 sigma). "Alcanzable" = 1 - sn^2 / MSE(const) con el campo
verdadero (el techo de cualquier término).
- **M2-a (positivo):** criterio: MSE drop >= 0.8 x alcanzable y cobertura en [0.90, 0.98].
- **M2-b (negativo):** sf = 0 (solo ruido). Criterio: |MSE drop| <= 0.03 y cobertura en [0.90, 0.98].
- **M2-c (descriptivo):** el mismo campo con ruido t de 3 grados (colas pesadas, como los picks de analista) de la
  misma varianza. Se reporta el cociente drop / alcanzable, para compararlo con el 0.6 medido en los datos reales.
Lectura: si M2-a y M2-b pasan, el código del GP da lo esperado y G-M2 se sostiene como fallo del gate en el piso de
ruido (como dice la nota del 14 de noviembre); si M2-a falla, G-M2 queda contaminado.

### C-BL. El blanqueo en el mecanismo de N-1, ACAP G-M5b, R0 y R1
Evento sin ruido hecho con `ttmt.radiation.radiation_operator` (7 estaciones al azar, tensor al azar), cada estación
dividida por una sd log-uniforme con cociente máximo/mínimo 7.6 (el medido en Ridgecrest). f = |Q^T v|^2 / |v|^2.
- **BL-a (negativo del defecto):** sd iguales, Q sin pesos. Criterio f >= 0.99.
- **BL-b (positivo del arreglo):** sd distintas, Q de diag(1/sd) G. Criterio f >= 0.99.
- **BL-c (descriptivo):** sd distintas, Q sin pesos (lo que hacen `N1-STACK/run.py` y `r0_run.Net`). Mediana sobre
  200 eventos, para comparar con el f sin ruido de `gcheck` de N-1 (0.705 Hi-net, 0.771 Ridgecrest).

### C-GEN. El desplazamiento fraccionario del generador (`ttmt.synth.inject._shift`)
Pulso de Ricker de 4 Hz a 100 Hz, retardo de 0.37 muestras. Criterio: error máximo frente al Ricker analítico
retardado <= 1e-3 del pico y energía dentro de 1e-3. Se reporta lo mismo con `np.interp` lineal (el patrón del
defecto de M6-SPLIT) como control negativo: debe perder energía de forma visible (> 1e-3).

## Revisiones estáticas (sin correr; cada una cita archivo y línea en la tabla)
- Qué guiones dividen por la sd de cada estación con un operador sin pesos (patrón 1).
- Qué guiones usan `np.interp` para un retardo (patrón 2).
- Qué tiempos de viaje usa cada línea (patrón 4).
- N-2: cómo se asocia cada forma de onda restada a su evento (`rebuild`: el orden de `comps` frente a `truth`).

## Clasificación
Cada línea queda como **se sostiene** (código comprobado con verdad conocida y sin defecto aplicable),
**contaminado** (defecto concreto, con archivo, línea y efecto medido o esperado) o **sin auditar todavía** (qué
falta y por qué: datos del portátil, arreglo sin fusionar). Las re-corridas con datos reales van después de fusionar
el arreglo del que dependen (#73, regla 21), con semillas nuevas, criterios escritos antes y sin leer ventanas
reservadas ni segmentos gastados (reglas 14, 16); se proponen al dueño en orden.

## Costo
Todo en CPU de la nube, minutos. No usa GPU, portátil ni datos del proyecto.

Resultados y clasificación: `tabla-auditoria.md` (24 de noviembre de 2026).
