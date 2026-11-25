# RAY1D: rayo en modelo 1-D en lugar de rayo recto (ficha, declarada antes de correr)

Diseño: supervisión, `supervision/tras-reunion-director-2026-11-20.md`, sección 2, copiado tal cual abajo. Lo que sigue
a la copia son los detalles de interfaz, fijados aquí antes de correr nada.

## Diseño (copia literal)

**Por qué.** El operador actual usa la línea recta fuente-estación para el patrón de radiación y para la polarización
en la estación. En un medio con velocidad que crece con la profundidad, el rayo sale de la fuente con otro ángulo y
llega a la superficie casi vertical, así que la P se ve sobre todo en la componente vertical. [hipótesis] Eso explica
parte de por qué el vector leído en eventos reales casi no se proyecta sobre el subespacio (f_emp 0.078 / 0.060).

**Operador G_ray.** Por nodo y estación, del modelo 1-D de cada conjunto (Ridgecrest: Hadley-Kanamori 1977, el
`PRED_MODELS` de `synth/realdata.py`; Hi-net: JMA2001 de `traveltime.py`):
- parámetro del rayo p de la tabla de tiempos (pendiente dT/dΔ);
- ángulo de salida en la fuente asin(p·v_fuente), con el acimut fuente-estación: dirección del rayo para el patrón
  de radiación P, SV y SH;
- ángulo de incidencia en la superficie asin(p·v_superficie): polarización P a lo largo del rayo que llega, SV y SH
  perpendiculares; las componentes S de la fuente se proyectan en SV y SH y se llevan a la estación;
- divergencia 1/d como hoy (no cambia el peso relativo entre estaciones más de lo que se mide aquí);
- pesado por 1/sd (ADR-0012).
Variante descriptiva: coeficientes de superficie libre de P y SV.

**Comprobación de código (nube, antes de todo).** Eventos fabricados con G_ray y leídos con G_ray, sin ruido: f ≥
0.99 en las geometrías reales de los dos conjuntos. Si no, no se corre nada más.

**Medidas y población (portátil).** Las de FE-FIX paso 3, en la banda operativa (0.5-6 Hz, 16 Hz, Ricker) y con Q
pesado: margen m con CAP6 en todos los objetivos con snr_med ≥ 5, y f_emp con g pesado; emparejadas G_ray frente a
G recto. Sin mitades: no se elige entre muchas opciones.

**Lectura (tres zonas, regla 36), sobre el cociente de medianas de m, G_ray / G recto, por conjunto:**
- ≥ +20 % en los dos conjuntos: G_ray entra al modelo (por ADR, antes del piloto);
- entre 0 y +20 %, o solo un conjunto sube: ganancia útil menor; decide el dueño;
- ≤ 0 en los dos: el rayo 1-D no basta; lo siguiente son las funciones de Green empíricas, como dijo el director.
Además, fijado ya: si f_emp sube al menos 0.2 en valor absoluto, la hipótesis de la polarización queda apoyada; si
sube menos de 0.05, no. Intervalos por realización, 10 000 remuestreos, semilla 20261125.

**Límites.** Rayo de alta frecuencia en un modelo 1-D; sin efectos de sitio ni dispersión; los EGF son la onda real
del evento original, así que ningún G sintético la representará del todo.

## Detalles de interfaz (fijados antes de correr)

1. **Código.** `src/ttmt/ray_operator.py` (`ray_operator(src, STA, model, free_surface=False)`), pruebas en
   `tests/test_ray_operator.py` (límite homogéneo igual a `radiation_operator`, P a lo largo del rayo que llega, S
   perpendicular). El modelo no lo usa salvo que un ADR lo adopte.
2. **Modelo y p.** p viene de `ttmt.synth.media.first_arrival` sobre `synth.realdata.pred_model(dataset)`, la misma
   llamada que predice los tiempos de llegada de los conjuntos. Para Hi-net es JMA2001 con los mismos puntos que
   `traveltime.JMA2001_POINTS`, en capas de 0.5 km en lugar de 0.25 km; la regla v2 prohíbe `ttmt.traveltime` en
   este camino. Rayo directo: sale hacia arriba; onda de cabeza: sale hacia abajo. Profundidad mínima 0.5 km, como
   en las predicciones.
3. **Comprobación de nube.** Etapa `check`: las 20 ubicaciones y los tensores de N1-STACK `gcheck2.json` (misma
   semilla y orden, verificados), geometría y modelo reales de cada conjunto, cadena operativa de N1-STACK (pre16,
   matched, pico de energía en el mejor desplazamiento común). La regla es la mediana de f ≥ 0.99 por conjunto en
   dos casos: sd = 1, y sd falsos (geomspace 1-10 permutado, semilla 20261126) con Q pesado. Descriptivo: los mismos
   sd falsos con Q sin pesar, la variante de superficie libre, y los eventos de G_ray leídos con el Q recto.
4. **Población y medida del portátil.** Etapa `measure`: todas las realizaciones de la partición tune, todos los
   objetivos con snr_med ≥ 5, sin mitades. Es el frente de referencia de FE-FIX paso 3: sd robusto del bloque de
   ruido, máscara de catálogo, llegadas predichas en la ubicación y el origen verdaderos.
   m = S_mt del objetivo / cuantil 0.999 del mismo S_mt sobre el ruido enmascarado sin desplazar, con el mismo G.
   CAP6 es el principal y sin tope es descriptivo. La superficie libre es descriptiva.
5. **f emparejada (pendiente de confirmar por el supervisor antes de correr).** f_emp de FE-FIX,
   (v·g_w)²/(‖v‖²‖g_w‖²), no depende de G, así que no puede cambiar entre G_ray y G recto. Se informa como ancla. La
   etapa calcula dos medidas que sí dependen de G, por conjunto, como diferencia de medianas G_ray − recto, con
   intervalo por realización:
   - f_g = ‖Q_wᵀ g_w‖²/‖g_w‖², la fracción de la onda empírica inyectada que cae en el subespacio. Sin ruido ni
     lectura, es la prueba directa de la polarización. **Medida de la regla 0.2 / 0.05 (confirmado por el supervisor; ver la corrección al final).**
   - f_an = ‖Q_wᵀ v‖²/‖v‖², el vector leído en el subespacio (descriptiva).
6. **Bootstrap.** Por realización (etiqueta), 10 000 remuestreos, semilla 20261125, una sola secuencia para los dos
   conjuntos en el orden de `--sets`. La ganancia es la mediana de m con G_ray / la mediana con G recto − 1.
7. **Veredicto.** El script imprime los números y las zonas de la ficha; la zona la lee el supervisor.

## Corrección del supervisor (20 de noviembre de 2026, antes de ver ningún valor de G_ray en datos)

Registrada también al final de `supervision/tras-reunion-director-2026-11-20.md`.
1. La regla 0.2 / 0.05 se aplica a f_g = ‖Q_wᵀ g_w‖²/‖g_w‖², leída como diferencia de medianas G_ray − recto, por
   conjunto. f_emp queda como ancla descriptiva y f_an como descriptiva. (f_emp no depende de G: error del diseño.)
2. Se acepta `pred_model` para el p de Hi-net. Se cuentan e informan por conjunto los rayos con p·v > 1 (seno
   recortado a 1), en la fuente o en la superficie.
Las zonas sobre m quedan como estaban escritas.

## Reglas 39-43 del dueño (22 de noviembre de 2026), aplicadas antes de correr la medida

- **39 y 40 (auditoría con verdad conocida y controles positivo y negativo).**
  - Control positivo de código: la etapa `check` (eventos de G_ray leídos con G_ray, f 0.999 / 0.998).
  - Control negativo de código: la misma etapa, con los eventos de G_ray leídos con el Q recto (0.925 / 0.866) y
    con los sd falsos sin pesar (0.716 / 0.797).
  - Controles de la medida, añadidos ahora a la etapa `measure` (semilla 20261127). Por cada objetivo, un evento
    sintético en su ubicación y su origen, hecho con G_ray (control positivo) o con el G recto (control negativo).
    Se suma al mismo ruido enmascarado y se escala para que su lectura sin ruido tenga el ‖v‖ del objetivo. Pasa
    por el mismo frente, la misma lectura, el mismo Q_w y el mismo cuantil de ruido.
    - El positivo pasa si f_g con G_ray es ≥ 0.99 y la ganancia de m de G_ray sobre el recto es > 0.
    - El negativo pasa si f_g con el recto es ≥ 0.99 y la ganancia es ≤ 0.
  - Si un control falla en un conjunto, la medida de ese conjunto no se lee y se informa al supervisor. En
    particular, un resultado negativo (zona ≤ 0) solo se declara si los dos controles pasan.
  - Las zonas y la regla de f_g no cambian.
  - Precisiones del supervisor (23 de noviembre, antes de correr):
    - (a) Cada criterio de control se aplica a la mediana por conjunto sobre los objetivos, el mismo estadístico
      que usan las zonas. Los controles se informan con su intervalo por realización, 10 000 remuestreos, semilla
      20261125, en un generador propio, así que la secuencia del bootstrap principal no cambia.
    - (b) El f_g ≥ 0.99 de cada control se calcula sobre el vector inyectado sin ruido. Es una comprobación del
      cableado y así se rotula.
    - (c) Si el control negativo falla solo porque la ganancia mediana de m queda por encima de 0, aunque sea
      poco y con un intervalo que incluya 0, se informa como fallado. Ese conjunto no se lee y no se reajusta nada.
- **41 (fusiones).** Después de cada fusión en main se trae main, se revisa la resolución y se declaran los archivos
  compartidos. Los de este PR, en el orden de fusión: `CHANGELOG.md`, `docs/artifact_manifest_and_error_log.md`
  (E66), `src/ttmt/radiation.py`, `src/ttmt/r1.py`.
- **42 (seguimiento).** `measure.json` y `measure_out.txt` se guardan en git desde el checkout del portátil, en la
  misma rama, antes de leer la zona.
- **43 (cierre).** Al cerrar, el informe dice qué reglas aplicaron y cómo.
