# N1-STACK: prueba previa de la suma con signo para eventos débiles (ficha)

16 de noviembre de 2026. Diseño del supervisor; lo ejecuta un hilo de implementación (Opus). Nada se ha corrido.
Origen: decisión del dueño en tarjeta, 16 de noviembre 11:39 AM Bogotá, "Priorizar eventos débiles", tras el
resultado de N-2 (PR #60).

Ficha subida al repositorio **antes** de correr (PR #63, tras la fusión del #60). El implementador fija detalles
de la interfaz del código; **no** cambia población, medida, umbral de ruido ni zonas sin volver al supervisor.

## Pregunta

De los eventos que el lector de fases no ve ni estando solos, ¿cuántos quedan por encima del ruido cuando se suma
la evidencia de toda la red con el signo que predice un mecanismo común, al mismo presupuesto de 30 falsas por hora?

Es el caso ideal (posición y tiempos verdaderos). Un número bajo cierra la idea antes de construir la búsqueda.
Un número alto solo dice que vale la pena construirla.

## Qué se sabe antes [medido]

- 120 objetivos T se pierden en el lector; con el vecino quitado, 88 siguen con menos de 6 llegadas (PR #60).
- El nivel de señal separa solo en parte recuperados de no recuperados (AUC 0.718 [0.593, 0.865]).
- Hi-net recupera 5 de 54 frente a 27 de 66 en Ridgecrest con el mismo nivel de señal. Causa sin probar.
- Techo teórico en C (ruido gaussiano, PR #47): un detector ideal de red recupera 0.75 de la pérdida de C.
  Es teoría con ruido ideal; esta prueba lo mide con el ruido real de los registros.

## Datos

Solo la partición de ajuste `b60d99b3d0bbbbfa` (los mismos 50 registros). Ninguna ventana reservada, ningún
segmento de evaluación gastado, ningún libro mayor. Ningún parámetro se ajusta con el resultado. Se anota como una
lectura exploratoria más de la partición de ajuste.

## Población

- **P-T (principal):** los 88 objetivos T × k que siguen con menos de 6 llegadas tras la resta exacta
  (`N2-REPICK/repick.json`, `n_repick < 6`). Se mide sobre el registro **con el vecino restado**, para aislar la
  poca señal del solapamiento.
- **P-C (secundaria):** los objetivos C con etapa "picker" en la línea base de PhaseNet+PyOcto (regla de STAGE-LOSS),
  sobre el registro sin modificar.
- **Control positivo:** objetivos que PhaseNet+PyOcto sí detecta (hasta 60, tomados en orden de etiqueta).
- Se reporta por conjunto (Ridgecrest, Hi-net) siempre, porque N-2 mostró que se comportan distinto.

## Medidas

Para cada objetivo, en su posición y tiempo de origen verdaderos, con los tiempos de llegada inyectados
(`arr_t_obs`) y una tolerancia de ±0.25 s en el tiempo de origen:

1. **S_mt (la nuestra):** el estadístico del detector de v2 (proyección sobre el subespacio de 6 componentes del
   tensor de momento, el mismo que usa el primer paso del modelo), sin pasar por el lector.
2. **S_en (referencia sin signo):** suma de energía por estación en las mismas ventanas, sin patrón de radiación
   (la idea de QuakeMigrate, un grado de libertad por estación). Sirve para saber si el signo aporta.

**Umbral, con ruido real.** Para cada registro, se recorre el bloque de ruido guardado (`noise_<bloque>`, sin
eventos inyectados) con el mismo estadístico, en todos los instantes y en una rejilla de n_loc = 100 posiciones
dentro del volumen de fuentes del conjunto. Se toman los máximos locales separados al menos 1 s. El umbral de cada
conjunto es el valor que esos máximos superan 30 veces por hora de ruido recorrido. Un umbral por conjunto y por
estadístico; se reporta cuántas horas de ruido lo sostienen. Sensibilidad: n_loc = 1000 en 5 registros.

**Número principal.** R = fracción de P-T con S_mt por encima del umbral. Intervalo: bootstrap sobre realizaciones
(10 000 remuestreos, semilla 20261122). Secundarios: lo mismo en P-C; lo mismo con S_en; diferencia S_mt − S_en
emparejada; todo por conjunto.

## Comprobaciones que detienen la corrida

- La línea base reproduce el objetivo congelado (T 0.2063, C 0.3115).
- Los eventos reconstruidos suman la señal guardada (la comprobación de N2-REPICK).
- Control positivo: al menos 0.90 de los objetivos que PyOcto detecta superan el umbral con S_mt. Si no, el
  estadístico o el umbral están mal y no se lee nada.
- El recorrido de ruido da entre 20 y 45 máximos por hora por encima del umbral al volver a contarlo (el umbral
  sale de esos mismos datos; es comprobación de aritmética, no validación).

## Lectura (regla 36: meta de ambición, tres zonas, sobre la estimación puntual de R en P-T)

| R | Lectura |
|---|---|
| ≥ 0.30 | Pasa la meta: se construye la fuente de candidatos por suma con signo (con búsqueda) |
| 0.10 a 0.30 | Meta no alcanzada, efecto útil: se reporta con intervalo y por conjunto; decide el dueño |
| < 0.10 | Baja: ni con posición y tiempos verdaderos la red ve esos eventos al presupuesto registrado; N-1 no se construye y se declara |

De dónde sale 0.30: la recuperación del mejor competidor en poca señal es 0.31 [medido]; recuperar 0.30 de lo que
pierde la sube en unos 0.20, por encima de la meta de +0.15 de la propuesta N-1. Los umbrales no se cambian después.

Lecturas adicionales fijadas de antemano, sin zona:
- Si S_mt − S_en no es positivo (intervalo que incluye 0), el signo no aporta en este régimen y se dice así: la
  ventaja sería de cualquier suma de red, no de nuestra idea.
- Si R en Hi-net es menor que la mitad de R en Ridgecrest, la causa de Hi-net no es solo el lector y hay que mirar
  las funciones de Green o el modelo de tiempos de ese conjunto antes de seguir.

## Límites

- Caso ideal: posición, tiempos y funciones de Green verdaderos. Es una cota superior. La búsqueda real paga más
  posiciones (umbral más alto) y errores de tiempo.
- Semisintético: funciones de Green empíricas y ruido real, no eventos débiles reales.
- El umbral usa pocas horas de ruido; 30 por hora queda estimado con pocas decenas de máximos. Se reporta el
  intervalo del umbral por remuestreo de registros.
- Superar el umbral no es una detección según la regla de acierto (posición y tiempo); no se afirma recuperación.

## Costo y reparto

- Horas de portátil (los datos de formas de onda están allí; entorno conda `ttmt-gpu`,
  `TTMT_DATA_ROOT=C:\Users\juanj\ttmt\realdata`, `TTMT_WORK_ROOT=C:\Users\juanj\ttmt\work`). Necesita la frase de
  autorización del dueño escrita en el hilo que lo corra, nombrando acción y carpeta.
- Implementación y corrida: Opus. Lectura del resultado y elección de rama: supervisor (Fable).
- Versión pequeña primero: 2 registros, sin lectura de zona.

## Tarea barata asociada: el lector en Hi-net

Con las salidas que ya existen (sin correr nada nuevo): comparar por conjunto, en los eventos **aislados** de la
partición, la fracción de llegadas verdaderas que PhaseNet lee frente al nivel de señal. Si a igual señal Hi-net lee
menos, la hipótesis del dominio de entrenamiento queda apoyada y el afinado del lector (N-3) debe incluir Hi-net.
Criterio: diferencia de fracción leída a señal igualada (tres tramos de señal), con intervalo por realización.

## Revisión del diseño antes de correr (16 de noviembre de 2026; nada se había corrido)

Motivo, hallado por el hilo de implementación en `data/manifests/synth`: la partición de ajuste tiene **un solo
bloque de ruido de 300 s por conjunto** (Ridgecrest 2100-2400 s, Hi-net 0-300 s; las 10 realizaciones de cada
conjunto lo comparten). Son 0.083 h por conjunto: a 30 por hora el umbral sería el 2.º o 3.er máximo local, y el
bloque contiene los eventos de catálogo (neutros en la verdad). Los demás bloques son de lectura o de EGF y no se usan.
Decisión del supervisor (mismo día, antes de cualquier corrida); población, medidas y zonas (0.30 / 0.10) no cambian:

1. Horas de ruido por **time-slides** sobre ese bloque: cada estación se desplaza circularmente 20 s × m_i, con los
   m_i distintos entre estaciones (diferencias de al menos 20 s), 100 deslizamientos, semilla fija. Antes de deslizar
   se **enmascaran por estación** las ventanas de llegada P y S de los eventos de catálogo del bloque; se reporta la
   fracción enmascarada.
2. Los time-slides rompen la coherencia del ruido entre estaciones (ya falló en M5: colas simultáneas 1.7 veces lo
   independiente), así que su umbral puede quedar bajo. **Umbral principal = el mayor de dos**: el de time-slides a
   30/h y el del bloque sin deslizar (su 3.er máximo local más alto, con la máscara; unos 30/h en 0.083 h). La zona
   se lee con ese máximo. Se reporta R con cada umbral por separado y el cociente de los dos umbrales; si el sin
   deslizar supera al de time-slides en más de 20 %, se dice como evidencia de ruido coherente.
3. El recuento de 20-45 máximos por hora se aplica solo al nulo de time-slides; es aritmética (da unos 30 por
   construcción), no una validación. El control positivo (≥ 0.90 de
   los aciertos de PyOcto sobre el umbral principal) se mantiene.
4. Límite: un solo bloque de 300 s por conjunto; el umbral sin deslizar descansa en 2-3 máximos y su intervalo es
   muy ancho.
5. El recorrido de ruido usa la misma regla de lectura que los objetivos, y la sd robusta por estación sale del
   bloque de ruido también para los registros con señal.

## Correcciones tras la revisión independiente de d5a93d1 (16 de noviembre de 2026, antes de correr)

1. La máscara de catálogo se aplica también al ruido de los objetivos: R = R_señal + R_ruido · (fuera de máscara)
   (− R_vecino en P-T); el t0 de cada registro es el del bloque (se comprueba y se detiene si no). Sin esto, el nulo
   pierde la energía de los eventos de catálogo y los objetivos no, y R sube. Sensibilidad: R con nulo y objetivos
   sin máscara.
2. Control positivo por conjunto: hasta 30 aciertos de PyOcto por conjunto (orden de etiqueta dentro del conjunto) y
   ≥ 0.90 en **cada** conjunto; si uno falla, no se lee nada.
3. Se reporta R también en los dos extremos del intervalo del umbral de time-slides.
4. La sd robusta por estación sale del bloque de ruido con sus eventos de catálogo; la infla algo, igual para nulo y
   objetivos, sin sesgo entre ellos.

## Primera corrida (NOCAP) y decisión del supervisor (16 de noviembre de 2026)

**Intento 1, registrado (commit 1950aeb; `baseline_out.txt`, `stats_smoke_out.txt`, `stats_out.txt`, `stats.json`).**
S_mt se programó sin tope (NOCAP), una elección de interfaz que se apartó del diseño ("el estadístico del detector de
v2"), que lleva el tope 6 registrado del MT-FW. La línea base se reprodujo (T 0.2063, C 0.3115). La corrida se
**detuvo en el control positivo**: Hi-net 0 de 6 y Ridgecrest 20 de 30, frente a 0.90 por conjunto. No se calculó
ni se leyó R. Lo que quedó a la vista antes de la decisión siguiente:
- Umbrales S_mt: Hi-net 76 748 con time-slides y 203 779 sin deslizar (cociente 2.66, fracción enmascarada 0.498,
  14 eventos de catálogo); Ridgecrest 222 y 324 (cociente 1.46, fracción enmascarada 0.298, 8 eventos).
- Medianas de S_mt: control Hi-net 893 y Ridgecrest 945; **P-T Hi-net 111 y Ridgecrest 40** (un vistazo a la
  población, aunque no se leyó R).

**Decisión del supervisor (mismo día, tras ver lo anterior; necesita el visto bueno del dueño por las reglas 14 y
20, pedido en tarjeta en el hilo del supervisor):**
1. Se repite con **CAP6** (cada bloque con |v| limitado a 6 antes de proyectar, `detector.CAP6`): vuelve al
   estadístico del diseño. S_en sigue sin tope. Población, medida R, umbral principal (el mayor de los dos) y
   zonas 0.30 / 0.10 no cambian.
2. Diagnóstico barato, con regla fijada antes: en los máximos de S_mt del nulo, la fracción de la energía leída
   (Σ_b |v_ib|², sin tope) que aporta la estación más fuerte, por conjunto: en el 2 % más alto de los máximos de
   cada deslizamiento y en los 10 máximos más altos sin deslizar. Si la mediana supera 0.5 **con NOCAP**, la
   hipótesis del transitorio de una sola estación queda apoyada; si no, la causa es otra y se dice. Se imprime
   también con CAP6.
3. Regla de parada con CAP6: control positivo **agrupado** ≥ 0.90 (como estaba escrito en el diseño); el control
   por conjunto se reporta, y un conjunto con menos de 10 (Hi-net tiene 6) se marca como insuficiente para
   juzgarlo solo. Esto reemplaza la regla por conjunto de las correcciones tras la revisión, y se declara como
   cambio tras ver datos. Si falla, la lectura es "instrumento no válido" y se para; no hay tercera variante sin
   el supervisor y el dueño.
4. Límites: en Hi-net la mitad del bloque queda enmascarada; el umbral sin deslizar es 2.66 veces el deslizado
   (ruido coherente o eventos sin catalogar). Se reporta el mismo cociente con CAP6.
5. R no se imprime hasta que se confirme la aprobación del dueño: `--read`.

Corridas: `stats ... --radial nocap --out stats_nocap_diag.json` (diagnóstico; repite el intento 1) y
`stats ... --radial cap6 --out stats_cap6.json` (sin `--read`).

**Aprobación del dueño** (tarjeta en el hilo del supervisor, 16 de noviembre de 2026, 4:48 PM Bogotá, opción "Sí,
repetir con el tope"): repetir N-1 con CAP6, mismo criterio y zonas, con la corrida sin tope registrada como intento.
Marca de la revisión independiente para cualquier lectura: CAP6 se fijó tras mirar la población (umbrales y medianas
de P-T del intento 1); su efecto en el control y en los umbrales se reporta junto a R.

## Resultado del intento 2 (CAP6) y del diagnóstico (16 de noviembre de 2026, portátil; commit 061d8d8)

Salidas: `stats_nocap_diag_out.txt`, `stats_nocap_diag.json` (NOCAP repetido, 440 s; mismas cifras que el intento 1)
y `stats_cap6_out.txt`, `stats_cap6.json` (495 s). 20 de 20 realizaciones; sin errores.

**Lectura por la regla declarada: con CAP6 el control positivo agrupado da 0 de 36 (Hi-net 0 de 6, Ridgecrest 0 de
30; hacía falta 0.90). Instrumento no válido; no se calculó ni se leyó R.** Por la decisión del supervisor, no hay
tercera variante sin él y sin el dueño.

| | Hi-net NOCAP | Hi-net CAP6 | Ridgecrest NOCAP | Ridgecrest CAP6 |
|---|---|---|---|---|
| umbral S_mt time-slides | 76 748 | 54.5 | 222 | 49.4 |
| umbral S_mt sin deslizar | 203 779 | 106.6 | 324 | 91.3 |
| cociente sin deslizar / time-slides | 2.66 | 1.96 | 1.46 | 1.85 |
| control sobre el umbral principal | 0 de 6 | 0 de 6 | 20 de 30 | 0 de 30 |
| mediana S_mt del control | 893 | 54.8 | 945 | 42.9 |

Diagnóstico de la estación más fuerte (regla: mediana > 0.5 con NOCAP):
- time-slides, 2 % más alto: Hi-net 0.992 y Ridgecrest 0.792. **Apoyada**: los máximos del nulo deslizado los
  pone una sola estación.
- 10 máximos más altos sin deslizar: Hi-net 0.304 (no apoyada) y Ridgecrest 0.528 (apoyada, por poco). En Hi-net
  el umbral sin deslizar, que es el principal, viene de energía en varias estaciones a la vez: ruido coherente o
  eventos sin catalogar.

Medida adicional, leída después de la corrida, sin regla (se reporta, no decide nada): la fracción de la energía
leída que cae en el subespacio del tensor, S_mt / S_en con NOCAP, en los aciertos de PyOcto es de mediana 0.259 en
Hi-net (6 eventos) y 0.207 en Ridgecrest (30). Con un vector al azar sería 6/42 = 0.14 y 6/66 = 0.09. El operador
de radiación analítico (rayos rectos, medio homogéneo, cocientes P/S fijos) explica así solo una quinta parte de los
eventos reales de EGF.
Hipótesis, sin probar: el estadístico falla sobre todo porque **G no se ajusta a las funciones de Green reales**, y
no solo por el ruido. Prueba que la cerraría: la misma fracción con un G empírico (los vectores g_ib medidos en
cada EGF, que inject guarda en `meta['targets'][...]['g']`), en los mismos aciertos.

**Lectura del supervisor (17 de noviembre de 2026):** N-1 queda como **instrumento no válido**. No se lee R ni hay
tercera variante. La fracción S_mt / S_en queda como dato descriptivo posterior, sin regla.

## Diagnóstico declarado, sin correr: G empírico frente a G analítico (17 de noviembre de 2026)

**Aprobado por el dueño** en tarjeta del supervisor (16 de noviembre de 2026, 9:27 PM Bogotá, "Sí, hacer el
diagnóstico"), solo en los aciertos de PyOcto. Declarado aquí antes de correr. No toca P-T ni R ni lee zona.

- Población: los mismos aciertos de PyOcto del control (`pops.json`, 30 Ridgecrest y 6 Hi-net), sobre el registro
  sin modificar, con la misma máscara, front end y lectura (pico en [-1, +8] muestras, desplazamiento común de −4..+4).
- G empírico: los vectores radiados de cada objetivo que guarda la inyección (`truth.meta['targets'][e]['g']`,
  forma (ns, 2, 3), Z cambiada a hacia abajo y ordenada como `radiation_operator`). Como un solo evento da un solo
  vector, el subespacio empírico tiene rango 1: f_emp = (v·ĝ)² / ‖v‖², con su referencia al azar 1 / (6 ns).
  G analítico: f_an = ‖Qᵀv‖² / ‖v‖² (rango 6), referencia 6 / (6 ns). Medidas emparejadas en cada objetivo, sin tope,
  por conjunto; se reporta también la diferencia emparejada.
- Lectura fijada antes de correr (supervisor), sobre la mediana de f_emp por conjunto: **> 0.6** la causa es G (la
  lectura del vector reproduce el evento, el modelo analítico no); **< 0.4** la causa está en la lectura del vector
  (ventana, pico, filtro) y se revisa antes que G; **0.4 a 0.6**, ambas.
- Límites: g se mide en velocidad filtrada 2-20 Hz a 100 Hz y v en la salida del filtro adaptado de 0.5-6 Hz a
  16 Hz; un cambio de signo por bloque entre las dos bandas baja f_emp, así que una f_emp baja no separa del todo
  "lectura" de "banda", y una f_emp entre 0.4 y 0.6 se lee como "no separa". El rango 1 frente a 6 hace que una f_emp alta sea la evidencia más fuerte.
- v: los vectores leídos en el desplazamiento común con mayor ‖v‖²; las dos fracciones sobre ese mismo v.
- Comprobación del código con un evento sintético analítico en ruido gaussiano (nube, sin datos del proyecto): f_an
  0.998 y f_emp 0.998 con g conocido, así que la lectura y el orden de filas reproducen un evento que sí sigue el
  modelo. Ese evento sintético no pasó por `inject` (pulsos Ricker con G analítico, Z hacia abajo invertida a Z
  hacia arriba como en los conjuntos): prueba la lectura y el orden de filas, no la convención de
  `meta['targets'][e]['g']`.
- Por eso (revisión de 2f33a2c), antes de leer la regla: g se pone a cero en los bloques sin llegada, igual que v; y
  se imprimen también f_emp sin girar Z y con P y S intercambiados. Si alguna de esas dos medianas supera a la de
  f_emp en más de 0.2 en un conjunto, la convención está mal y la zona de ese conjunto no se lee.
- Uso: `run.py gdiag pops.json --out gdiag.json` (programado, sin correr).

## Control de la cadena con G analítico (supervisor, 17 de noviembre de 2026; declarado antes de correr)

Dentro del diagnóstico que aprobó el dueño (no toca P-T ni lee zona). Propuesta de la revisión independiente,
aprobada por el supervisor.

- 10 eventos por conjunto con el **G analítico** (`radiation_operator` × tensor al azar, semilla 20261125), en las
  posiciones verdaderas de los aciertos de PyOcto del control, con llegadas del modelo de predicción del conjunto.
  Pulso Ricker de 4 Hz en velocidad en cada llegada, Z hacia arriba en los datos, sumado al bloque de ruido del
  conjunto como lo hace `inject` (`_taper` de 0.5 s y desplazamiento fraccionario `_shift`), en un origen al azar
  del bloque. Escala: la que da a S_en la mediana de S_en de los aciertos del conjunto (`stats_nocap_diag.json`).
- Mismo front end, máscara y lectura; f = S_mt / S_en sin tope, sobre la lectura con mayor ‖v‖².
- Lectura fijada: mediana de f **≥ 0.8** la cadena (convenciones, lectura, filtro) está bien y una f baja en los EGF
  es física; **< 0.5** hay un error de convención o de lectura y se arregla antes de interpretar nada; entre ambas,
  se reporta y no se concluye.
- Límite: no usa `inject.realise` completo (que planea sus propios pares y escalas por SNR de tarjeta), sino sus
  mismas primitivas de suma; el pulso es analítico, no un EGF.
- Comprobación en la nube con un lector falso (ruido gaussiano, 7 estaciones): mediana 0.969. Uso:
  `run.py gcheck pops.json stats_nocap_diag.json --out gcheck.json`, en la misma sesión que `gdiag`.

## Resultado del diagnóstico y del control de la cadena (17 de noviembre de 2026, portátil; commit 028aa3f)

Salidas: `gdiag_out.txt`, `gdiag.json`, `gcheck_out.txt`, `gcheck.json`. Lectura por las reglas declaradas arriba;
la interpretación final es del supervisor.

| | Hi-net (n 6, < 10) | Ridgecrest (n 30) |
|---|---|---|
| f_emp mediana [q25, q75] (azar) | 0.291 [0.258, 0.369] (0.024) | 0.033 [0.020, 0.127] (0.015) |
| f_an mediana (azar) | 0.135 (0.143) | 0.129 (0.091) |
| f_emp − f_an emparejado, mediana | +0.159 | −0.095 |
| regla gdiag (> 0.6 G, < 0.4 lectura) | < 0.4: lectura | < 0.4: lectura |
| convención: sin girar Z / P-S intercambiados | 0.316 / 0.000: coherente | 0.037 / 0.001: coherente |
| gcheck, G analítico, S_mt/S_en mediana [mín, máx] | 0.685 [0.245, 0.819] | 0.768 [0.707, 0.838] |
| regla gcheck (≥ 0.8 bien, < 0.5 error) | entre ambas: no se concluye | entre ambas: no se concluye |

Lo que se mide, sin ir más allá de las reglas:
- Por la regla de `gdiag`, la causa señalada es **la lectura del vector**, en los dos conjuntos: el vector leído
  apenas se parece al g propio del evento (Ridgecrest 0.033, casi el azar 0.015). Por el límite declarado, esa
  regla no separa "lectura" de "banda" (g a 2-20 Hz, v a 0.5-6 Hz).
- `gcheck` queda en la zona media: con el G analítico en los registros reales, la cadena recupera 0.69 y 0.77 de la
  energía en el subespacio, frente a 0.97 con ruido gaussiano en la nube. No se concluye que la cadena esté bien
  ni que tenga un error de convención; la comprobación de convención de `gdiag` tampoco señala error.
- f_an de los aciertos es aquí 0.13 (mismo v para las dos fracciones), menor que el 0.21-0.26 de S_mt/S_en de
  `stats` (allí cada estadístico toma su propio máximo sobre desplazamientos).

## Comprobación de banda (supervisor, 17 de noviembre de 2026; declarada antes de correr)

Aprobada por el dueño en la tarjeta del supervisor el 17 de noviembre de 2026, 7:16 AM Bogotá ("Sí, comprobar la banda"); frase de Remote Control escrita en el hilo a las 7:19 AM. No toca P-T ni R.

- Población y g: los mismos aciertos y el mismo g de `gdiag` (bloques sin llegada en cero, Z hacia abajo).
- Front end elegido: **el mismo filtro con que se midió g**, `realdata.bandpass` 2-20 Hz (Butterworth de fase cero)
  a 100 Hz, sin filtro adaptado ni remuestreo; cada estación dividida por la sd robusta de su bloque de ruido
  filtrado igual; máscara de catálogo llevada a 100 Hz. Misma regla de lectura en segundos: pico de energía en
  [−0.06, +0.5] s de cada llegada, desplazamientos comunes de −0.25..+0.25 s en pasos de 1/16 s, la lectura con
  mayor ‖v‖².
- Medida: f_emp en 2-20 Hz, emparejada con f_emp en 0.5-6 Hz (front end de R1, recalculada en la misma corrida).
- Lectura fijada: mediana por conjunto **> 0.6** en 2-20 Hz, la causa es la banda; **< 0.4**, la causa es la regla
  de lectura (pico o ventana) y no la banda; entre ambas no separa.
- Descriptivo, sin regla: razón de S_en (normalizado por el ruido de cada banda) entre 2-20 Hz y 0.5-6 Hz por acierto.
- Comprobación del código en la nube (sin datos del proyecto): un evento analítico con pulso de 8 Hz y g medido
  por `inject.radiated_vectors` (el camino real de los conjuntos) da f_emp 0.998 en 2-20 Hz y 0.566 en 0.5-6 Hz.
  Esto confirma además la convención de `meta['targets'][e]['g']` que usan `gdiag` y `gband`.
- Uso: `run.py gband pops.json --out gband.json`.

### Comprobaciones añadidas tras la revisión de a7cea2c (declaradas antes de correr, descriptivas, sin regla)

La revisión independiente señaló dos dudas sobre los resultados de `gdiag`/`gcheck`: en Ridgecrest f_an (0.129)
supera a f_emp (0.033), lo que podría indicar un desajuste de índices entre v y g; y el 0.69-0.77 de `gcheck` no
separa la cadena del ruido. Se anotan como límites de esos resultados y se añaden tres medidas, que no cambian
ninguna regla ni ningún umbral:

- `gdiag` y `gband` paran si el orden de estaciones de la verdad (`meta['stations']`) difiere del manifiesto.
- `gdiag` y `gband` imprimen la correlación de rangos (Spearman) entre las normas por bloque |v| y |g| sobre los
  bloques con g no nulo. Cerca de cero o negativa haría sospechar un desajuste de índices. En la nube: +0.996
  (`gband`, evento sintético) y +0.837 (`gdiag`).
- `gcheck` imprime f sobre la lectura sin ruido (f_clean). En la nube: 0.998 (f con ruido 0.969). Si en el portátil
  f_clean es alto y f con ruido queda en 0.69-0.77, la pérdida es del ruido y no de la cadena.
- Se corren junto con `gband`, si el dueño lo aprueba, repitiendo `gdiag` y `gcheck` con las mismas semillas.

## Resultado de la comprobación de banda y de las comprobaciones de la revisión (17 de noviembre de 2026, portátil; commit 2508ddf)

Corridos con las mismas semillas sobre e435c6a, sin cambiar código ni reglas. El control de orden de estaciones se
hizo de verdad en los 36 aciertos (gband y gdiag) y coincide con el manifiesto. `gcheck` se corrió con el ref
documentado (`stats_nocap_diag.json`), igual que la vez anterior. Salidas: `gband*`, `gdiag2*`, `gcheck2*`
(`gdiag.json` y `gcheck.json` quedan intactos).

| Medida | Hi-net (n 6) | Ridgecrest (n 30) | Regla declarada |
|---|---|---|---|
| f_emp 2-20 Hz, mediana [q25, q75] | 0.569 [0.481, 0.664] | 0.597 [0.480, 0.706] | entre 0.4 y 0.6: no separa |
| f_emp 0.5-6 Hz (misma corrida) | 0.291 | 0.033 | — |
| diferencia emparejada, mediana | +0.208 | +0.509 | descriptivo |
| razón de S_en 2-20 / 0.5-6 Hz | 4.26 | 2.62 | descriptivo |
| correlación de rangos \|v\|-\|g\| en 2-20 Hz | +0.720 | +0.670 | descriptivo |
| gdiag repetido: f_emp / f_an | 0.291 / 0.135 | 0.033 / 0.129 | igual que 028aa3f: "lectura del vector" |
| gdiag: correlación de rangos \|v\|-\|g\| en 0.5-6 Hz | +0.714 | +0.421 | descriptivo |
| gcheck repetido: S_mt/S_en con ruido | 0.685 | 0.768 | entre 0.5 y 0.8: sin conclusión |
| gcheck: f sobre la lectura sin ruido | 0.705 | 0.771 | descriptivo |

Lectura según la ficha, sin veredicto propio: la comprobación de banda **no separa** en ninguno de los dos
conjuntos (Ridgecrest 0.597 queda justo bajo 0.6; Hi-net tiene n < 10). La lectura de la zona le toca al supervisor.
Observaciones descriptivas, sin regla previa: las correlaciones de rangos son positivas en ambas bandas, lo que no
apoya un desajuste de índices entre v y g; y en `gcheck` f sin ruido es casi igual a f con ruido, de modo que la
pérdida del 0.23-0.30 no viene del ruido sino de la propia cadena con G analítico (en la nube, con lector falso,
f sin ruido fue 0.998).

## Lectura del supervisor y cierre de N-1 (17 de noviembre de 2026)

- Por la regla declarada, `gband` queda en la zona media en los dos conjuntos (0.569 y 0.597): **no separa**.
- Medido, descriptivo emparejado: leer en 2-20 Hz sube f_emp en +0.21 (Hi-net) y +0.51 (Ridgecrest), y la energía
  normalizada en 2-20 Hz es 2.6 a 4.3 veces la de 0.5-6 Hz. Con `gcheck` sin ruido en 0.705 y 0.771 hay además una
  pérdida de unos 0.25 en la propia cadena de lectura con registros reales, que el lector falso de la nube no tenía.
- Lectura: **dos causas a la vez**, la banda (efecto grande) y una pérdida de la cadena de lectura; ninguna sola
  explica todo.
- **Hipótesis abierta, sin probar:** la explicación de la pérdida de la cadena (sin ruido, con G analítico) no se ha
  medido.
- N-1 se cierra aquí; no se corre nada más. Dónde se arregla el front end lo decide el dueño en el hilo del supervisor.

## Implementación (`run.py`; detalles de interfaz fijados antes de correr y aceptados por el supervisor)

- Front end de R1 (`datasets._preprocess`): media y tendencia fuera, Butterworth 4 de fase cero 0.5-6 Hz, 16 Hz,
  Z hacia abajo; filtro adaptado Ricker 2 Hz de `ttmt.kernel`; cada estación dividida por la sd robusta de su salida
  en el bloque de ruido. Todo es lineal, así que señal, ruido y vecino se filtran por separado y se suman.
- Lectura: en cada bloque (estación, fase), el vector de 3 componentes en el pico de energía dentro de [-1, +8]
  muestras de la llegada (−0.06 a +0.5 s, como `inject.radiated_vectors`); llegada ausente = vector cero.
  Tolerancia de origen: máximo sobre desplazamientos comunes de −4..+4 muestras, cada estadístico el suyo.
- S_mt = ½‖Qᵀv‖² en el intento 1; desde la decisión de arriba, ½‖Qᵀφ(v)‖² con φ el tope 6. Q base ortonormal de `radiation_operator` en la posición verdadera con todas las estaciones del
  registro (`ProjectedEnergy` sin tope). S_en = ½‖v‖² en las mismas lecturas.
- Recorrido de ruido: 100 posiciones uniformes en la caja de las fuentes objetivo del conjunto (profundidad ≥ 0.5 km,
  semilla 20261124), tiempos del modelo de predicción del conjunto (`realdata.pred_model`: HK1977 / JMA2001);
  máximo sobre posiciones en cada muestra de origen; máximos locales separados ≥ 1 s (supresión de no máximos), lo
  que contiene el máximo sobre ±0.25 s. Máscara de catálogo: [tP − 1 s, tS + 8 s] por estación, tiempos del mismo
  modelo. Umbral de time-slides: el valor superado floor(30 × horas) veces; intervalo por remuestreo de
  deslizamientos (2000). Sensibilidad n_loc = 1000 frente a 100 en los 5 primeros deslizamientos (el bloque es uno
  solo, así que "5 registros" se lee como 5 deslizamientos).
- P-T: filas de `N2-REPICK/repick.json` con `n_repick < 6`, sobre el registro con el vecino restado (reconstrucción
  de N2-REPICK, que se detiene si los eventos no suman la señal guardada). P-C y control: etapa `baseline`
  (reproduce el objetivo congelado o se detiene; control: 30 por conjunto). Bootstrap de R por realización, 10 000, semilla 20261122.
- Versión pequeña: `stats --limit 2` (una realización de cada conjunto, nulo completo), solo comprobaciones, sin R
  ni lectura.
- Prueba con datos falsos en la nube (ruido gaussiano, 7 estaciones, 20 deslizamientos): recuento 29.4 por hora;
  sin señal S_mt 12, con señal creciente sube como el cuadrado de la amplitud. No toca ningún dato del proyecto.

Orden en el portátil (entorno `ttmt-gpu`, `TTMT_DATA_ROOT`, `TTMT_WORK_ROOT` como N2-REPICK):

    python experiments/exploratory/N1-STACK/run.py baseline C:/Users/juanj/ttmt/tune/ss_ridgecrest_2019_egf_v1/tune_phasenet.json C:/Users/juanj/ttmt/tune/ss_hinet_2012_egf_v1/tune_phasenet.json --out pops.json
    python experiments/exploratory/N1-STACK/run.py stats pops.json experiments/exploratory/N2-REPICK/repick.json --limit 2 --out stats_smoke.json
    python experiments/exploratory/N1-STACK/run.py stats pops.json experiments/exploratory/N2-REPICK/repick.json --out stats.json

La zona la lee el supervisor; este hilo solo entrega las salidas.

## Nota posterior (FE-FIX, 17 de noviembre de 2026; descriptivo)

La revisión de FE-FIX mostró que `gdiag` comparaba v/sd con g sin pesar. Recalculado con g y Q pesados por 1/sd
(`FE-FIX/run.py gdiagw`, commit 7f1cc44): f_emp 0.078 en Hi-net (n 6) y 0.060 en Ridgecrest (n 30), frente a 0.291 y
0.033 sin pesos.

**Corrección (supervisor, 18 de noviembre de 2026):** la explicación de banda de la lectura del supervisor de N-1 era
incompleta. Con g pesado, f_emp baja en lugar de subir, así que el blanqueo no explica la lectura baja. Hipótesis sin
probar: el G analítico de rayos rectos no representa las señales reales (EGF con efectos de sitio y dispersión).
