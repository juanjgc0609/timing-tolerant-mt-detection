# AUDIT-NEG: tabla de auditoría de los negativos (regla 39), 24 de noviembre de 2026

Criterios y controles declarados antes de correr en `card.md` (commit 3a2d667). Salida: `controls_out.txt`
(nube, CPU, numpy 2.5.3, scipy 1.18.1, 1083 s; semilla 20261124). Los 13 controles con criterio pasan, así que
los instrumentos de la auditoría son válidos (regla 40). Ningún número registrado cambia. Ningún criterio se relaja.

## Controles

| control | resultado | criterio | lectura |
|---|---|---|---|
| GEN-a `inject._shift` | error 1.0e-14, energía +0 | <= 1e-3 | PASA |
| GEN-b `np.interp` lineal (control negativo) | energía −1.8e-2 | > 1e-3 | PASA (detecta el patrón de M6-SPLIT) |
| BL-a sd iguales, Q sin pesos | f mín 1.000 | >= 0.99 | PASA |
| BL-b sd distintas, Q de diag(1/sd) G | f mín 1.000 | >= 0.99 | PASA |
| BL-c sd distintas (7.6), Q sin pesos | f mediana 0.799 [0.732, 0.857] | descriptivo | mismo orden que el gcheck de N-1 sin ruido (0.705 / 0.771) |
| M2-a campo conocido + ruido gaussiano | caída 0.401 de 0.423 alcanzable (0.95), cobertura 0.953 | >= 0.8 x alcanzable, 0.90-0.98 | PASA |
| M2-b solo ruido | caída −0.005, cobertura 0.950 | abs <= 0.03, 0.90-0.98 | PASA |
| M2-c campo + ruido t3 | caída 0.335 de 0.356 (0.94), cobertura 0.954 | descriptivo | el GP llega al 94 % de lo alcanzable incluso con colas pesadas |
| M5-a fondo exacto, H 200 000 | producto 4.3e-4 / 3e-5; S_avg 2.1e-4 / 2e-5 (p 1e-2 / 1e-3); media exp(S_avg) 0.979 | <= p; [0.9, 1.1] | PASA |
| M5-b fondo demasiado ligero | producto 5.3e-2 a 1e-3 | > 1e-3 | PASA |
| M5-c potencia | P(S_avg >= log 1000) 0.755 | >= 0.5 | PASA |
| M5-d lecturas de pico, fondo radial con n_fit 30 000 y 2 000 000 | S_avg y producto 0 a 1e-3 en los dos | descriptivo, lectura fijada | ver M5-EVALUE abajo |

## Clasificación por línea

| línea | negativo registrado | patrón que aplica (archivo:línea) | evidencia | clasificación |
|---|---|---|---|---|
| **N-1** (N1-STACK) | "instrumento no válido": control positivo CAP6 0 de 36 | blanqueo: `N1-STACK/run.py:394-395` divide por la sd de cada estación; `:130-132` y `:402` proyectan sobre `basis(radiation_operator)` sin pesos | gcheck sin ruido 0.705 / 0.771, las mismas cifras del defecto medido por FE-FIX (PR #68); BL-c lo reproduce con verdad conocida (0.80) y BL-b muestra que el arreglo lo lleva a 1.000 | **contaminado (blanqueo)**. También la banda 0.5-6 Hz (gband, ya registrado). R nunca se leyó: la población P-T está intacta |
| **R0** (Ridgecrest, registrado) | −0.383 [−0.481, −0.310] | blanqueo: `code/r0/r0_run.py:50-51` y `:65-67` (Gs, Qst sin pesos); además la búsqueda de ±1 s de los términos (D1 / E50) y el medio homogéneo VP 6.0 (`code/r0/r0_lib.py:5`) | defectos en el código; efecto en el número sin medir. La saturación del tope (C18c, AUC 0.39-0.49) es independiente del blanqueo | **contaminado (blanqueo, D1, medio homogéneo)**. El registrado sigue registrado; la hora EVAL está gastada, así que solo cabe una re-corrida descriptiva. Hipótesis sin probar: el tamaño (−0.38) hace poco probable que cambie de signo |
| **R1** (Hi-net, registrado) | ver `results/taskR1` | blanqueo: `src/ttmt/r1.py:60-61` y `:71`. Tiempos JMA2001 1-D (`r1.py:33`), no de rayo recto | defecto en el código, efecto sin medir | **contaminado (blanqueo)**. Igual que R0: descriptivo sobre datos gastados |
| **C18a** | "ningún medio predice mejor las lecturas" | búsqueda de ±1 s truncada (E50). El blanqueo no aplica: lee tiempos, y escalar una estación no mueve su pico | E50, ya registrado | **contaminado (E50); la conclusión ya estaba retirada** |
| **M5-ACAP** | ACAP retirado: G-M5a falla en las dos semillas; G-M5b falla en EVAL (N4 0.616) | G-M5a usa la cadena 2M (`code/t2m_lib.py:328`, sd escalar por ventana): sin blanqueo. G-M5b usa `r0_run.Net`: blanqueo | G-M5a: MT@S0 reproduce las 72 de 72 y 49 de 49 partes registradas (control positivo de la cadena); NOCAP queda último en el día con transitorios (control negativo) | **se sostiene** por G-M5a, que basta para la decisión. G-M5b queda **contaminado (blanqueo)** como lectura secundaria; EVAL está gastado y no se rehace |
| **M5-EVALUE, pasos 1-2** (producto) | falla 1e-3 en ruido real (5.8e-3 / 1.2e-3) | sin blanqueo (sd escalar, `run_c.py:104-108`). El código del e-valor pasa con verdad conocida (M5-a, b, c), también a 1e-3, que los tests no cubrían | ninguna variante del fondo arregla el exceso; el diagnóstico señala colas simultáneas entre estaciones (1.7 veces lo independiente) | **se sostiene**: el código es correcto y el fallo viene de los datos (dependencia entre estaciones) |
| **M5-EVALUE, paso 3** (S_avg) | 1.4e-3 a 1e-3 en ruido falso i.i.d. (cota 95 % 2.4e-3; 28 800 hipótesis) | `--local` (blanqueo) solo en una variante secundaria | M5-d: con lecturas de pico ideales, el fondo radial da S_avg = 0 a 1e-3, con n_fit 30 000 y con 2 000 000. La lectura declarada ("error finito del fondo" o "defecto en el código") no aplica: **el exceso no se reproduce fuera de la cadena** | **sin auditar todavía**: `m5.py` queda absuelto, pero la causa está en la cadena de `run_c.py` (lector sobre la salida correlacionada del filtro adaptado, hipótesis de red que comparten picos, ajuste por estación) y no se ha localizado. Además el exceso es marginal (la cota inferior del intervalo no está calculada) |
| **M2-GP** (G-M2) | caída del MSE 0.191 < 0.30 en confirmación | tiempos HK 1-D con rayos directos y refractados (`src/ttmt/traveltime.py:59`), **no** rayo recto homogéneo; el GP modela además el residuo | M2-a y M2-b pasan: el código recupera el 95 % de lo alcanzable y no inventa estructura en ruido puro | **se sostiene** como fallo del gate en el piso de ruido (nota del 14 de noviembre). Dato nuevo, descriptivo: en sintético el GP llega al 94-95 % de lo alcanzable; en los datos reales, al ~60 % de su propia cota, así que la cota de los datos reales está sobreestimada o el residuo real es menos suave que un Matern 3/2 |
| **N-2** (N2-REPICK) | neto −8 [−23, 7], zona baja | revisado `rebuild` (`N2-REPICK/run.py:100-137`): las formas restadas se guardan en el orden de los objetivos, y en `inject.realise` los objetivos van primero en la verdad (`src/ttmt/synth/inject.py:259-272` y la construcción de `Truth`), así que `comps[f]` es el evento f. El desplazamiento es el de frecuencia (GEN-a) | controles que pasaron en la corrida: la suma reconstruida igual a la guardada (<= 1e-3), el picker reproduce los picks guardados, K_local = K en 120 de 120 | **se sostiene en lo revisado**. Falta un control con verdad conocida: restar el propio objetivo debe dejarlo sin picks. Es barato y necesita el portátil. Defecto del control K por debajo de 0.5 s ya registrado en la ficha |
| **M1A-ONSET** | no es un negativo (G-M1a pasa en EVAL) | usa `r1.Net` (blanqueo), pero solo lee tiempos de pico y de PhaseNet, que no cambian al escalar una estación | revisión estática | **no aplica** |

## Patrones buscados en todo el código

- Blanqueo (división por la sd de cada estación con un Q sin pesos): `code/r0/r0_run.py`, `src/ttmt/r1.py`,
  `code/explore_contrib/c18a_timing.py` (solo tiempos), `N1-STACK/run.py`, `M5-EVALUE/run_c.py --local`. Las
  simulaciones M6-* usan `node_operators` con ruido gaussiano de igual sd (BL-a: sin efecto si las sd son iguales).
- Retardo por `np.interp`: además del generador de M6-SPLIT, `experiments/exploratory/M6-THEORY/confounders.py:164`.
  Es de la línea M6, fuera de este alcance: se avisa a su hilo. El generador de los conjuntos semisintéticos
  (`inject._shift`) es exacto.
- Rayo recto homogéneo: R0 (VP 6.0) y el operador de amplitudes `radiation_operator` (hipótesis de FE-FIX, que prueba
  RAY1D en #73). M2 y R1 usan tiempos 1-D.
- Búsqueda local con un solo arranque (patrón de `free_one_event`): no aparece en estas líneas. `GP1.fit` ya arranca
  desde tres escalas.

## Re-corridas propuestas, en orden (el dueño decide; ninguna corre antes de fusionar #73)

1. **N-1 con el arreglo del blanqueo** (después de #73; portátil). Mismo control positivo agrupado >= 0.90, mismas
   zonas 0.30 / 0.10, población y partición de ajuste intactas. Es el negativo con más probabilidad de cambiar:
   el defecto explica la pérdida de la cadena que se midió. La banda se decide aparte (ADR-0012 la deja igual).
2. **N-2, control "restar el propio objetivo"** (portátil, ~3 min, GPU). Cierra la auditoría de N-2. No hace falta #73.
3. **M5 paso 3: localizar el exceso en la cadena** (nube, ruido falso). Mismo `run_c.py avg --weighted` con el fondo
   ajustado sobre una muestra independiente grande del mismo ruido falso. Si el exceso sigue, la causa es la cadena
   (dependencia entre hipótesis o entre P y S) y se busca; si desaparece, es el error finito del ajuste. Criterio a
   escribir antes de correr.
4. **R0, R1 y ACAP G-M5b, descriptivos con el arreglo** (portátil, sobre datos gastados). No cambian ningún registro.
   Solo dicen si "negativo" sigue siendo la palabra correcta para la tesis. Prioridad baja.

No hace falta re-correr: M5-ACAP (G-M5a), M5-EVALUE pasos 1-2, M2-GP, C18a (ya retirado) y M1A.
