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
     lectura, es la prueba directa de la polarización. **Propuesta como medida de la regla 0.2 / 0.05.**
   - f_an = ‖Q_wᵀ v‖²/‖v‖², el vector leído en el subespacio (descriptiva).
6. **Bootstrap.** Por realización (etiqueta), 10 000 remuestreos, semilla 20261125, una sola secuencia para los dos
   conjuntos en el orden de `--sets`. La ganancia es la mediana de m con G_ray / la mediana con G recto − 1.
7. **Veredicto.** El script imprime los números y las zonas de la ficha; la zona la lee el supervisor.
