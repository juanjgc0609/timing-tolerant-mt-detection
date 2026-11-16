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

## Implementación (`run.py`; detalles de interfaz fijados antes de correr y aceptados por el supervisor)

- Front end de R1 (`datasets._preprocess`): media y tendencia fuera, Butterworth 4 de fase cero 0.5-6 Hz, 16 Hz,
  Z hacia abajo; filtro adaptado Ricker 2 Hz de `ttmt.kernel`; cada estación dividida por la sd robusta de su salida
  en el bloque de ruido. Todo es lineal, así que señal, ruido y vecino se filtran por separado y se suman.
- Lectura: en cada bloque (estación, fase), el vector de 3 componentes en el pico de energía dentro de [-1, +8]
  muestras de la llegada (−0.06 a +0.5 s, como `inject.radiated_vectors`); llegada ausente = vector cero.
  Tolerancia de origen: máximo sobre desplazamientos comunes de −4..+4 muestras, cada estadístico el suyo.
- S_mt = ½‖Qᵀv‖², Q base ortonormal de `radiation_operator` en la posición verdadera con todas las estaciones del
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
