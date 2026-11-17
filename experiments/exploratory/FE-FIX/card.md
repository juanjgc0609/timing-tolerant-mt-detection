# FE-FIX: arreglar la lectura del modelo antes del piloto (diseño)

17 de noviembre de 2026. Diseño del supervisor; lo implementa un hilo de programación (Opus); la lectura de cada
paso la hace el supervisor. Origen: decisión del dueño en tarjeta, 17 de noviembre 7:45 AM Bogotá, "Arreglar antes
del piloto", tras N-1 (PR #63).

Ficha de FE-FIX (copia del diseño del supervisor, `/mnt/project-files/supervision/fe-fix-lectura-del-modelo-2026-11-17.md`).
Se sube **antes** de correr cada paso. El implementador fija la interfaz; no cambia población, medidas ni umbrales sin volver al supervisor.

## Qué se sabe [medido, PR #63]

- La lectura actual (front end R1: 0.5-6 Hz, remuestreo a 16 Hz, Ricker de 2 Hz) da f_emp 0.291 (Hi-net) y 0.033
  (Ridgecrest) en los aciertos de PyOcto; a 2-20 Hz sin Ricker, 0.569 y 0.597. La energía en 2-20 Hz es 2.6 a 4.3
  veces la de 0.5-6 Hz.
- Con eventos fabricados con el G analítico y **sin ruido**, la cadena real da 0.705 y 0.771; el lector simplificado
  de la nube daba 0.998. Hay una pérdida de unos 0.25 dentro de la cadena, de causa desconocida.
- Candidata propuesta por la revisión: solape P/S en estaciones cercanas (S-P menor que la ventana de lectura).

## Datos

Solo la partición de ajuste `b60d99b3d0bbbbfa`. Nada reservado, ningún segmento gastado, ningún libro mayor. Lo que
se elige aquí (banda, ventana) se elige en una mitad y se confirma en la otra.

## Paso 1. Solape P/S (descriptivo, sin correr nada nuevo)

Sobre `gcheck2.json`: f sin ruido frente al S-P mínimo del evento; y f sin ruido quitando los bloques con S-P < 0.6 s.
Lectura fijada: si al quitarlos la mediana sube a 0.90 o más, el solape P/S es la causa principal; si sube menos de
0.05, no lo es; entre medias, es parte.

## Paso 2. Dónde se pierde la cadena (sin ruido, eventos analíticos de gcheck)

Partir del lector de la nube (f 0.998) e ir cambiando **un elemento a la vez** hacia la cadena real, midiendo f sin
ruido en cada paso: geometría y estaciones reales; remuestreo a 16 Hz; filtro 0.5-6 Hz; Ricker de 2 Hz; giro y signo
de componentes; ventana de lectura; regla del pico y desplazamiento común. Se reporta la tabla completa.
Lectura: la causa es el elemento que más baja f; si ninguno baja más de 0.05 por sí solo, la pérdida es de
interacción y se dice. Meta de la cadena arreglada: f sin ruido ≥ 0.95 en los dos conjuntos.

## Paso 3. Elegir la banda y confirmarla

- Población: todos los objetivos de la partición (T, C e aislados) con `snr_med` ≥ 5 de la verdad inyectada, no solo
  los aciertos de PyOcto (Hi-net tenía 6). Se parten por realización en mitad de elección y mitad de confirmación,
  con semilla fija 20261123, estratificado por conjunto.
- Candidatas, fijadas ya: 0.5-6 Hz (actual), 1-10 Hz, 2-20 Hz, 2-30 Hz; cada una con la cadena arreglada del paso 2.
- Se elige la de mayor mediana de f_emp en la mitad de elección (la misma banda para los dos conjuntos; si difieren
  por más de 0.1, se reporta y decide el dueño).
- Lectura en la mitad de confirmación, sobre la mediana de f_emp por conjunto (regla 36, tres zonas):
  ≥ 0.6 pasa; 0.4 a 0.6 meta no alcanzada, efecto útil, decide el dueño; < 0.4 baja: el G analítico no basta y hay
  que volver a la física (funciones de Green empíricas o corrección por estación) antes del piloto.
- Intervalo: bootstrap por realización, 10 000 remuestreos, semilla 20261124.

## Después

- El cambio de front end es un cambio de diseño del modelo: se registra en ADR-0012 (reservado por el plan v2.18) con estos
  números antes del piloto.
- La prueba N-1 queda cerrada; repetirla con la lectura nueva sería una ficha nueva que decide el dueño.

## Límites

- Semisintético con funciones de Green empíricas; el G empírico tiene rango 1 y se midió a 2-20 Hz, lo que favorece a
  esa banda: por eso la elección se confirma también con el gcheck analítico sin ruido, que no depende de g.
- Elegir y confirmar en la misma partición de ajuste, en mitades distintas; no es evidencia confirmatoria.

## Costo

Paso 1, minutos en la nube. Pasos 2 y 3, horas de portátil (frase de autorización del dueño en el hilo que corra).

## Paso 1: implementación y resultado (17 de noviembre de 2026, nube, sin datos de forma de onda)

- `run.py sp ../N1-STACK/gcheck2.json --out sp.json`. Para cada evento analítico de `gcheck2.json` calcula el S-P de
  cada estación con el modelo de predicción del conjunto en la posición verdadera (la misma llamada `first_arrival`
  que usó `gcheck`; estaciones del manifiesto en git). Solo usa geometría: no lee ningún registro.
- `gcheck2.json` no guarda la lectura por bloque. Si algún evento tuviera bloques con S-P < 0.6 s, quitarlos exigiría
  una corrida nueva y el script se detiene. Si no hay ninguno, f no cambia por construcción.
- Descriptivo, añadido por el implementador y sin regla: bloques con S-P < 0.75 s (fin de la ventana de lectura,
  +8 muestras a 16 Hz, más el mayor desplazamiento común, +4) y Spearman entre el S-P mínimo y f sin ruido.
- Orden de los hechos: la regla se declaró en el diseño del supervisor (`fe-fix-lectura-del-modelo-2026-11-17.md`,
  17 de noviembre de 2026, 7:47 AM Bogotá) antes del cálculo; el cálculo (geometría determinista sobre una salida que
  ya existía) se hizo en la nube y esta sección se subió después.

| Conjunto | Eventos | S-P mínimo por evento, mediana [mín, máx] | Bloques con S-P < 0.6 s | f sin ruido, antes → sin ellos | Bloques con S-P < 0.75 s | Spearman (S-P mín, f) |
|---|---|---|---|---|---|---|
| Hi-net | 10 | 1.91 s [1.77, 1.92] | 0 de 70 | 0.705 → 0.705 (+0.000) | 0 | +0.33 |
| Ridgecrest | 10 | 2.08 s [1.44, 2.68] | 0 de 110 | 0.771 → 0.771 (+0.000) | 0 | −0.53 |

Por la regla de la ficha (subida < 0.05), **el solape P/S no es la causa** en ninguno de los dos conjuntos. Ninguna
estación está a menos de 1.4 s de S-P, y la ventana de lectura termina a 0.75 s como máximo. La lectura le toca al
supervisor. Salidas: `sp.json`, `sp_out.txt`.

## Paso 2: implementación (17 de noviembre de 2026; declarada antes de la corrida en el portátil)

`run.py chain ../N1-STACK/gcheck2.json [--laptop] --out chain.json`. Mismos 10 eventos analíticos por conjunto que
`gcheck` (tensores y orígenes resorteados con su semilla 20261125 y su orden; el script se detiene si un origen no
coincide con `gcheck2.json`), sin ruido. Cada configuración da la mediana [mín, máx] de f = ‖Qᵀv‖²/‖v‖².

**Corrección de una premisa del diseño.** El 0.998 de la nube no salió de un lector simplificado: era la misma
función `gcheck` (mismo remuestreo, filtro, Ricker, giro, ventana, pico y desplazamiento) con un lector falso de
datos: 7 estaciones uniformes en ±20 km a cota 0, modelo de Ridgecrest y ruido gaussiano de la misma sd en todas las
estaciones. Por eso la tabla tiene dos partes:

1. **Hacia la cadena real**, un elemento a la vez desde la comprobación de la nube: geometría y estaciones del
   conjunto (manifiesto en git y posiciones de `gcheck2.json`); modelo de velocidad del conjunto; normalización por
   estación con la sd robusta del bloque de ruido del conjunto (necesita el portátil). También la versión acumulada.
2. **Desde la cadena real**, quitando un elemento a la vez, con la alternativa fijada aquí: sin normalización por
   estación (sd = 1); sin remuestreo (lectura a 100 Hz, misma ventana y desplazamientos en segundos, Ricker de 2 Hz
   muestreado a 100 Hz); sin el pasabanda 0.5-6 Hz; sin el filtro adaptado; sin el giro de Z (control de signo, no
   debe subir f); ventana más ancha [−0.25, +1.0] s; regla de polarización principal de la ventana (vector propio
   mayor del segundo momento 3×3, signo del pico de energía, norma √λ) en lugar del pico; sin desplazamiento común.

Lectura fijada por el supervisor: la causa es el elemento que más baja f; si ninguno baja más de 0.05 por sí solo,
es interacción. Meta de la cadena arreglada: f sin ruido ≥ 0.95 en los dos conjuntos.

**Prueba del código en la nube (17 de noviembre, antes de subir esta sección; solo geometría en git, ningún
registro):** con las estaciones, posiciones y modelo reales y sd = 1, f sin ruido es 0.999 [0.997, 0.999] en
Hi-net y 0.998 [0.996, 0.999] en Ridgecrest, igual que la comprobación de la nube. Todo lo que no depende del ruido
del conjunto queda en ~0.998; el único elemento que falta medir para llegar a 0.705 / 0.771 es la normalización por
estación. Salida: `chain_cloud_out.txt`.

**Hipótesis (sin probar hasta la corrida en el portátil):** dividir cada estación por su sd saca a v del espacio de
G, porque la base Q se calcula con G sin ponderar. Si se confirma, el arreglo natural es ponderar también el
operador (Q de diag(1/sd)·G), y la pregunta de si el modelo principal tiene el mismo desajuste es del supervisor.

**Comprobación declarada por el supervisor (17 de noviembre, antes de correr):** (a) f sin ruido con la
normalización por sd_i y Q de G sin pesos (la cadena actual) y (b) lo mismo con Q de diag(1/sd)·G (configuración
`qweight`). Lectura: si (b) sube a 0.95 o más y (a) se queda en 0.70-0.77, la causa es el blanqueo y el arreglo es
(b); si (b) no sube más de 0.05 sobre (a), no lo es. Se reporta el cociente máximo/mínimo de sd_i por conjunto. Si se
confirma, el arreglo entra en ADR-0012 junto con la banda. El supervisor señala que `r1.Net` (src/ttmt/r1.py:60-61)
divide por `station_scale` y construye Gs y Qst sin pesos: el modelo principal tendría la misma estructura (hipótesis).

Prueba del código de (b) en la nube con sd falsas (lognormales, semilla 5, máximo/mínimo 11.7) y geometría de Hi-net:
(a) 0.847, (b) 1.000. Solo prueba que el código hace lo que dice; no dice nada de las sd reales.

Detalle de interfaz (corregido tras un fallo de dimensiones en la primera corrida, 4770cb8, que se detuvo en Ridgecrest antes
de normalizar y no se guardó): en la configuración con geometría falsa (7 estaciones) y sd del conjunto se usan las sd de las
7 primeras estaciones del manifiesto, en orden. En Hi-net (7 estaciones) no cambia nada.

Costo en el portátil: menos de 1 minuto, CPU, sin GPU; lee solo el bloque de ruido de la partición de ajuste.

## Paso 2: resultado (17 de noviembre de 2026, portátil; commit 24b2e4a, 31 s, CPU)

| f sin ruido, mediana [mín, máx] | Hi-net | Ridgecrest |
|---|---|---|
| comprobación de la nube | 0.999 [0.997, 0.999] | 0.999 [0.997, 0.999] |
| + geometría del conjunto / + modelo del conjunto | 0.999 / 0.999 | 0.998 / 0.999 |
| + normalización por estación (sola) | 0.837 [0.753, 0.910] | 0.834 [0.658, 0.939] |
| (a) cadena real (gcheck) | 0.705 [0.624, 0.819] | 0.771 [0.728, 0.849] |
| (b) cadena real con Q de diag(1/sd)·G | **0.999 [0.998, 1.000]** | **0.999 [0.997, 1.000]** |
| real sin normalización (sd = 1) | 0.999 | 0.998 |
| real sin remuestreo / sin pasabanda / sin Ricker | 0.707 / 0.705 / 0.690 | 0.779 / 0.771 / 0.770 |
| real sin giro de Z (control de signo) | 0.196 | 0.435 |
| real con ventana ancha / polarización principal / sin desplazamiento | 0.705 / 0.707 / 0.705 | 0.771 / 0.778 / 0.771 |
| cociente máximo/mínimo de sd_i | 5.70 | 7.63 |

Por las reglas declaradas: el elemento que más baja f es la normalización por estación (−0.29 y −0.23), y (b) sube
a 0.999 en los dos conjuntos con (a) en 0.70-0.77. **Por la regla del supervisor, la causa es el blanqueo y el arreglo
es (b)**, que cumple la meta de la cadena (≥ 0.95). Ningún otro elemento mueve f más de 0.02, salvo el control de
signo, que la baja como debe. La lectura le toca al supervisor.

**Lectura del supervisor (17 de noviembre de 2026):** por la regla declarada, la causa es el blanqueo y el arreglo es
(b): 0.999 / 0.999 frente a 0.705 / 0.771. Ningún otro elemento mueve f más de 0.02, y el control de Z invertida baja
como debe.

## Paso 3: implementación (17 de noviembre de 2026; declarada antes de correr)

`run.py band --weight --out band.json`. Lo que fija el implementador, sin tocar población, medidas ni umbrales:

- Población: todos los objetivos (`truth.target`) de todas las realizaciones de ajuste de cada conjunto con
  snr_med ≥ 5, donde snr_med es la mediana sobre estaciones de `truth.snr` (SNR por estación de la verdad inyectada).
  Incluye T, C y aislados; el compañero de un par queda en el registro.
- Mitades: por realización y por conjunto, permutación con semilla `[20261123, suma de los códigos del id del conjunto]`;
  la primera mitad (redondeo hacia arriba) elige, la otra confirma.
- Front end de cada banda candidata (0.5-6, 1-10, 2-20 y 2-30 Hz): quitar media y tendencia, Butterworth 4 de fase cero
  a 100 Hz, Z hacia abajo; **sin remuestreo ni filtro adaptado**, porque 16 Hz no puede contener bandas por encima de
  8 Hz y el Ricker de 2 Hz pertenece a la banda de 0.5-6 Hz (el paso 2 mostró que ninguno de los dos mueve f). Cada
  estación se divide por la sd robusta de su bloque de ruido en esa banda; máscara de catálogo de N1-STACK llevada a 100 Hz.
- Lectura: la regla de N1-STACK en segundos (pico de energía en [−0.06, +0.5] s, desplazamiento común ±0.25 s).
- Medida: f_emp = (v·g_w)²/(‖v‖²‖g_w‖²) con el g inyectado del propio objetivo (Z hacia abajo, bloques ausentes en
  cero) y **g_w = g/sd por estación** (el arreglo del paso 2, pendiente de la lectura del supervisor; sin `--weight`
  se usa g). Descriptivos: f_an con Q de diag(1/sd)·G, y la cadena actual R1 (16 Hz, Ricker) como referencia.
- Elección: la banda con mayor media de las dos medianas por conjunto en la mitad de elección; si la mejor banda de
  un conjunto supera a la elegida en más de 0.1 en ese conjunto, se reporta y decide el dueño.
- Confirmación: mediana de f_emp por conjunto en la otra mitad con la banda elegida; tres zonas (≥ 0.6 pasa; 0.4-0.6
  meta no alcanzada, efecto útil, decide el dueño; < 0.4 baja). Intervalo 95 % por bootstrap de realizaciones,
  10 000 remuestreos, semilla 20261124.
- Límite (de la ficha): f sin ruido de los eventos analíticos de `gcheck` en cada banda, con el arreglo, independiente de g.
- Prueba del código en la nube con un lector falso (geometría de Hi-net, sd distintas por estación, pulsos de 8 Hz,
  SNR alto): todas las bandas dan 1.000, la cadena R1 0.996-0.999. Solo prueba el código.
- Costo previsto en el portátil: unos minutos, CPU, sin GPU; lee solo la partición de ajuste.
- Confirmado por el supervisor (17 de noviembre) con dos condiciones: `--weight` (g y Q pesados por 1/sd) y, como
  referencia descriptiva que no entra en la elección, la banda 0.5-6 Hz con la cadena operativa del modelo (16 Hz y
  Ricker) con el mismo peso (fila «R1 operative chain»).
- Límite: la frecuencia de muestreo y el filtro operativos del modelo los fija ADR-0012 después, con su costo de
  cómputo. El Ricker se eligió por la potencia de detección en ruido, que f no mide.

### Revisión del paso 3 antes de correr (supervisor, 17 de noviembre de 2026, 10:34 AM Bogotá)

Motivo: la revisión independiente de d39221a señaló que elegir la banda por f_emp favorece 2-20 Hz por
construcción, porque g se midió en 2-20 Hz. El control analítico tampoco sirve para elegir: sin ruido da ~1.0 en
todas las bandas, y con ruido su pulso fabricado decidiría la banda. Es un cambio de diseño anterior a cualquier
corrida; lo decide el supervisor y no necesita al dueño. **Sustituye la medida de elección y la confirmación de la
sección anterior.** Población, mitades, bandas, front end, lectura y semillas no cambian.

- **Medida de elección: margen de detección m** de cada objetivo. m = S_mt del objetivo (posición verdadera, llegadas
  observadas, la lectura de N1-STACK con desplazamiento común) / cuantil 0.999 del mismo S_mt sobre todos los orígenes
  del bloque de ruido **sin deslizar**, con los eventos de catálogo enmascarados, en la misma posición (llegadas del
  modelo de predicción del conjunto, orígenes cada 1/16 s). S_mt = ½‖Q_wᵀφ(v)‖² con Q_w de diag(1/sd)·G (el arreglo
  del paso 2). Fijado por el implementador: φ = CAP6, el estadístico registrado; sin tope como descriptivo.
- Elección: la banda con mayor media de las medianas de m por conjunto en la mitad de elección. Si la mejor banda de
  un conjunto supera en ese conjunto a la elegida en más del 10 %, se reporta y decide el dueño.
- **Confirmación** en la otra mitad, emparejada contra la referencia operativa (0.5-6 Hz, 16 Hz, Ricker, con pesos).
  Ganancia = mediana de m con la banda elegida / mediana de m de la referencia − 1, por conjunto. Tres zonas: **+20 %
  o más** pasa y se adopta; **entre 0 y +20 %** la ganancia es útil pero menor y decide el dueño; **0 o menos** se
  queda la banda actual (con pesos) y se dice. Intervalo 95 % por bootstrap de realizaciones, 10 000 remuestreos,
  semilla 20261124.
- f_emp con g/sd queda solo como chequeo de lectura de la banda elegida: por debajo de 0.4 se avisa y no se adopta sin
  mirar. **f_emp está sesgado hacia 2-20 Hz por construcción** (g medido en 2-20 Hz) y no dice qué banda detecta mejor.
- Ítem 4 (descriptivo): `run.py gdiagw ../N1-STACK/pops.json --out gdiagw.json` recalcula el f_emp de `gdiag` de N-1 con
  g y Q pesados por 1/sd; el resultado se anota también en la ficha de N-1.
- Ítem 5: el supervisor verificó leyendo el código que `r1.Net` divide por `station_scale` (r1.py:60-61), que
  `node_operators` usa solo geometría y que `p1_windows` divide además MF por su sd. La prueba ejecutada va en las
  pruebas unitarias de ADR-0012.
- Prueba del código en la nube (lector falso, geometría de Hi-net, sd distintas, pulsos de 8 Hz, SNR alto): `band`
  corre de punta a punta. Con SNR alto el tope CAP6 satura todos los bloques, y los márgenes quedan entre 4.2 y 5.9.
  `gdiagw` da f_emp pesado 0.997. Solo prueba el código.
- Uso: `run.py band --out band.json` y `run.py gdiagw ../N1-STACK/pops.json --out gdiagw.json`. Costo previsto: unos
  minutos, CPU, sin GPU.
