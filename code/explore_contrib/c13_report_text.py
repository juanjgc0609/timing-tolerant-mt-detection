"""Prose of REPORT_C13.txt (Spanish, ASCII).  PARTS: strings and table references resolved by c13_report.py."""

P_HEAD = r"""
C13 - CONTROL DE FALSOS DESCUBRIMIENTOS POR BLANCO-SENUELO (TARGET-DECOY / NULA EMPIRICA) PARA EL DETECTOR MT-FW
(explore_contrib) - INFORME
======================================================================================================================
Estado: EXPLORATORIO, NO CONFIRMATORIO. Nada de lo que sigue prueba una hipotesis registrada.
Fecha: lunes 5 de octubre de 2026, 17:03-18:50 (-05). Maquina de 2 nucleos; un solo proceso de calculo a la vez
(OPENBLAS_NUM_THREADS=1, OPENBLAS_CORETYPE=Prescott), mas llamadas cortas de prueba con prioridad baja.
Tiempo de pared: lectura y codigo 17:03-17:08; corrida principal 17:08:36-18:02:32 (54 min); senuelo corregido D3c
18:15:50-18:26:04 (10 min); gemelos sin eventos 18:26-18:35; analisis (bootstrap) 18:35-18:41; ademas 13 min de una
primera version de gemelos y analisis (18:02-18:15) que se descarto (seccion 2); redaccion el resto.
No se edito el prerregistro ni ningun documento de /mnt/user-data/outputs; no se modifico ningun archivo anterior de
explore_contrib (solo archivos nuevos c13_*).

Semilla y datos (SOLO semilla base 999; ninguna otra semilla fue generada, leida ni tocada):
  redes     task2a_confirm.network(999, celda, j) con j = 1400-1411 (12 redes por celda y regimen), celdas ladder_4
            (4 estaciones), ladder_6 (6), unilateral_c40 (8, cobertura de un solo lado); regimenes thin y strong;
            estaticas sd 0.3 s; correcciones +SC del catalogo registrado de 20 eventos de cada red.
            Indices usados antes en explore_contrib: 0-4, 10-15, 20-35, 40-89, 99, 130-149, 160-179, 195, 199-239.
            grep previo de 14xx como indice de red o TAG en *.py, *.sh, *.txt: sin apariciones.
  TAGs      flujos [999, TAG, ci, j, ri, 41|233|61], todos nuevos (usados antes: 2003, 3001-3008, 3701-3706,
            3801-3812, 3901-3913, 4000-4004, 4101-4105):
            1401 ventanas STD (xc.draw_sources, 0-2 eventos, kmax 5), 20 por red;
            1402 ventanas DENSE (c2_run.draw_dense, 3-5 eventos en 9 km, kmax 7), 20 por red;
            1403 ventanas CLUT (sin eventos, clutter del regimen, kmax 5), 8 por red;
            1404 ventanas NOISE (sin eventos y sin clutter, kmax 5), 2 por red;
            1405 sorteos de los senuelos ([999, 1405, ci, j, ri, tipo de ventana, ventana], y [..., 7] para D3c);
            1406 prueba de humo (solo red 195, la red de validacion de C10/C11).
  Muestra   72 partes (2 regimenes x 3 celdas x 12 redes), todas terminadas, ninguna excluida: 1440 ventanas STD
            (1178 eventos; 40% de las ventanas sin evento), 1440 DENSE (5747 eventos), 576 CLUT, 144 NOISE.
            El generador gen_realistic.py no se toco; se usa por c8_lib.generate_site, como en C9-C11.
Detector: MT-FW con tope por estacion 6 (c10_lib.mt_loop, modo hard, +-0.25 s, correcciones +SC registradas, reticula
de 3 km x 3 profundidades = 675 nodos, origenes 5.0-8.0 s paso 0.1, exclusion 5 km, deflacion voraz), sin cambios. Es
"MT@S0" de C11. Puntuacion de un atomo = estadistico de seleccion; registro = minimo acumulado (como en el FROC exacto).
Reglas de acierto: LOC+T (epicentro a <= 5 km Y |error de origen| <= 0.5 s) y LOC (epicentro a <= 5 km).
"""

P_DESIGN = r"""
----------------------------------------------------------------------------------------------------------------------
1. PLAN DEL EXPERIMENTO (codigo: c13_lib.py, c13_run.py, c13_run_d3c.py, c13_twin.py; analisis c13_analyze.py)
----------------------------------------------------------------------------------------------------------------------
Idea. El detector se ejecuta dos veces sobre cada ventana: una con las hipotesis reales (BLANCO) y otra con hipotesis
SENUELO que no pueden corresponder a una fuente real. Si las detecciones falsas del blanco se comportan como las del
senuelo, el numero de senuelos por encima de un umbral estima el numero de falsos del blanco por encima de ese umbral,
sin etiquetas. Cada ventana recibe un sorteo nuevo de cada senuelo (flujo propio, TAG 1405). Senuelos:
  D1p  moveout con estaciones permutadas: la tabla de tiempos de la estacion i es la de la estacion pi(i), pi un
       desarreglo al azar (ningun punto fijo); el termino de estacion sigue a la estacion fisica; operadores de
       amplitud (G, Q) sin cambio, de modo que la estructura de Gram por nodo es la misma.
  D1o  moveout con desfases: cada columna (estacion, fase) de la tabla de muestras se desplaza con su propio desfase
       circular sobre los 40 s de traza, desfases separados al menos 1 s entre si; G, Q sin cambio.
  D2s  plantillas con signo invertido: la mitad de los bloques (estacion, fase) de 3 componentes cambia de signo en G
       y Q; tiempos sin cambio.
  D2r  plantillas con polarizacion revuelta: una rotacion 3D al azar por estacion aplicada a los bloques de ambas fases
       en G y Q; tiempos sin cambio. (La inversion temporal de la plantilla no es un senuelo aqui: el detector lee
       trazas filtradas con una ondicula simetrica.)
  D3   sustituto de datos: cada traza de estacion se desplaza circularmente con su propio retardo al azar sobre toda
       la traza (retardos separados al menos 3 s); hipotesis sin cambio. Equivale a D1o con un desfase por estacion
       (conserva el par P-S de cada estacion).
  D3c  D3 con retardos restringidos a [-7.5 s, +3.5 s] y separados al menos 1 s (POST HOC, ver seccion 2 y 4).
Calibraciones del umbral sobre un lote de ventanas (blancos y senuelos agrupados en el lote):
  tdc  busqueda separada: se aceptan los blancos con registro >= t*, t* el menor t con
       (1 + #senuelos >= t) / max(1, #blancos >= t) <= q.
  bh   p-valor conformal p = (1 + #senuelos >= s) / (n_senuelos + 1) de cada atomo blanco y Benjamini-Hochberg al
       nivel q sobre los m atomos blancos del lote. Como n_senuelos = m, es algebraicamente tdc con el factor
       m / (m + 1); se calculo y da los mismos numeros que tdc en todas las lineas (tabla TH).
  comp competencia pareada por rango (knockoff+): atomo r del blanco contra atomo r del senuelo de la misma ventana.
  fpw  (agregado durante la corrida, seccion 2) estimacion por senuelos de los eventos falsos por ventana: umbral =
       puntuacion de senuelo superada por a lo sumo c x (ventanas del lote) senuelos; es la version desplegable del
       punto de operacion del benchmark (<= c eventos falsos por ventana).
  Niveles: q = 0.05, 0.1, 0.2 (y 0.3, 0.5 como referencia); c = 0.1, 0.25.
Lotes: 20 ventanas = las 20 ventanas de un tipo de una red (72 lotes por tipo de ventana); 120 ventanas = 6 redes de
la misma celda y regimen (particiones al azar de las 12 redes en dos grupos; estimacion puntual promediada sobre 20
particiones). Limitacion: el lote de 120 mezcla 6 geometrias de red de la misma celda, no son 120 ventanas de una red.
Metricas: FDP de un lote = atomos aceptados falsos / max(1, aceptados); FDR realizada = media de la FDP sobre lotes
(un lote que no acepta nada cuenta 0); eventos falsos por ventana realizados; recall = eventos distintos acertados /
eventos; "oraculo" = el mejor umbral unico elegido con las etiquetas sobre las mismas ventanas con a lo sumo el mismo
numero de falsos (para fpw, tambien el punto del benchmark a <= c). Un duplicado de un evento ya acertado no cuenta
como falso ni como acierto nuevo (como en el FROC exacto).
Incertidumbre: bootstrap jerarquico (redes dentro de celda x regimen, luego ventanas dentro de red; 400 replicas;
umbrales, lotes y punto del oraculo recalculados en cada replica), intervalos percentiles del 95%. Los p-valores de
intercambiabilidad se obtienen permutando, ventana por ventana, las listas de blanco y de senuelo (1000 permutaciones).
"""

P_DECL = r"""
----------------------------------------------------------------------------------------------------------------------
2. QUE SE FIJO ANTES DE VER RESULTADOS Y QUE ES POST HOC
----------------------------------------------------------------------------------------------------------------------
Fijado antes de lanzar la corrida principal (17:08:36; sha256 de c13_lib.py fa0cb11a...267eef y de c13_run.py
56695e84...4149cf en ese momento): el detector, los cinco senuelos D1p, D1o, D2s, D2r, D3 y sus constantes, las reglas
tdc, bh y comp, los niveles q, las dos reglas de acierto, los dos tamanos de lote, redes y TAGs. No se cambio ninguna
constante de detector, senuelo o calibracion despues de ver datos; no se ajusto nada para que un metodo gane.
Miradas intermedias (para probar el codigo de analisis): con la red 1400 (17:12) y con las redes 1400-1403 (17:28).
POST HOC o agregado despues de esas miradas:
  (a) fpw: agregada a las 17:15, con una red vista, antes de tener la muestra. Es la pregunta de los auditores puesta
      en la misma moneda que el benchmark (eventos falsos por ventana), no una variante elegida por su resultado.
  (b) Gemelos sin eventos (c13_twin.py): planeados con una red vista. La PRIMERA version estaba mal: copiaba los
      generadores aleatorios con copy.deepcopy, que pierde la secuencia de semillas y su contador de hijos, y el
      generador de formas de onda engendra un hijo por ventana; el "gemelo" contenia restos de los eventos. Se detecto a
      las 18:18 al comparar las puntuaciones maximas del gemelo con las de ventanas vacias de la misma red, se corrigio
      (clon con secuencia de semillas; cada ventana se regenera completa desde clones y se compara muestra a muestra) y
      se repitio todo. Las salidas invalidas se conservan como out/c13_twin_parts_BUGGY_v1 y no entran en ningun numero.
  (c) D3c: agregado a las 18:15 despues de ver que D3 y D1o NO son intercambiables con el blanco en ventanas de clutter
      fuerte SIN eventos (seccion 4). Causa identificada en el generador: el clutter se sortea entre 2 s y 36 s de la
      traza y las hipotesis leen las muestras 95-451 (5.9-28.2 s); un retardo circular arbitrario trae los extremos sin
      clutter a la zona leida y el senuelo ve menos clutter. D3c restringe los retardos para que no ocurra. D1o tiene el
      mismo defecto y NO se repitio (falta de tiempo): sus numeros en clutter fuerte llevan ese sesgo.
  (d) Tabla TE2 (intervalos media +- 1.96 errores estandar entre lotes independientes): agregada al ver que los
      intervalos percentiles del bootstrap de la FDR quedan desplazados hacia arriba respecto de la estimacion puntual
      (remuestrear ventanas con reemplazo reduce la diversidad efectiva del lote y un umbral elegido con los datos
      aprovecha mas el azar). Las dos formas se dan; la marca * usa el intervalo que indica cada tabla.
  (e) Lineas agrupadas por celda ("both ladder_4", etc.): agregadas al ver que la celda unilateral queda
      sistematicamente por encima de las otras.
Sin correccion por multiplicidad: 6 senuelos x 4 calibraciones x 2 lotes x 5 niveles x 2 reglas x 18+6 lineas. Todo es
descriptivo.
No hecho: D1o corregido; BOMP-PK; el refinamiento T3 y los terminos JHD de C11 (el detector es MT@S0 duro); lotes de
120 ventanas de una misma red; calibracion conformal con un conjunto de calibracion separado; senuelos sobre el
residuo deflactado del blanco; ventanas con tasas de eventos mas bajas que STD.
Comprobaciones del codigo: (1) los motores senuelo D2s y D2r conservan Q ortonormal (desvio 1e-15) y el rango de G;
(2) la tabla de muestras reconstruida en c13_lib coincide con la del motor de C10 (assert por red); (3) c13_twin.py y
c13_run_d3c.py regeneran las ventanas y reproducen exactamente los atomos blanco guardados (assert en cada parte);
(4) en ventanas de solo ruido los seis senuelos dan la misma distribucion que el blanco (seccion 4).
"""

P_SUM = r"""
----------------------------------------------------------------------------------------------------------------------
0. RESUMEN Y VEREDICTO
----------------------------------------------------------------------------------------------------------------------
(1) NINGUN senuelo da control valido en TODAS las celdas con las dos reglas de acierto y con potencia util. El
    resultado depende de la regla de acierto y de la construccion del senuelo.
(2) En ventanas SIN eventos todos los senuelos bien construidos son intercambiables con el blanco (D1p, D2s, D2r, D3c:
    p de permutacion entre 0.04 y 1.0 en los 6 estratos, uno solo de 48 por debajo de 0.05; cociente de colas
    0.75-1.33). D3 y D1o, tal como se fijaron, NO lo son en clutter fuerte (D3: p 0.001-0.011 y 1.5-1.75 veces mas
    blancos que senuelos en la cola): defecto del retardo circular sobre toda la traza, corregido en D3c.
(3) En ventanas CON eventos los falsos del blanco no son "clutter": los generan los propios eventos. Frente al gemelo
    de la misma ventana sin eventos (mismo clutter y ruido), el numero de atomos falsos por ventana por encima del
    umbral t95 es 12 veces mayor en STD [9, 18] y 45 veces mayor en DENSE [33, 63]; a mas de 10 km de cualquier evento,
    6.7 [4.9, 10.3] y 24 [17, 37] veces. Por eso un senuelo solo sirve si conserva la energia de los eventos dentro de
    las muestras leidas y reproduce coincidencias parciales entre estaciones.
(4) Regla LOC: D1p (estaciones permutadas) controla. tdc con D1p, FDR realizada en las lineas agrupadas (lote 120 /
    lote 20): q = 0.05 -> 0.02-0.05 / 0.00-0.02; q = 0.1 -> 0.07-0.10 / 0.02-0.08; q = 0.2 -> 0.14-0.16 / 0.13-0.16.
    Ningun estrato individual viola q = 0.2. Reserva: la celda unilateral_c40 queda por encima del nominal a q <= 0.1
    con lotes de 120 (DENSE: 0.10 [0.07, 0.16] a q = 0.05 y 0.15 [0.11, 0.19] a q = 0.1); con lotes de 20, 0.04 y 0.12.
    D2r tambien controla bajo LOC con algo menos de recall; D2s controla pero no acepta casi nada; D3c queda en 1.3-2
    veces el nominal (0.10 / 0.16-0.18 / 0.21-0.26); D3 y D1o fallan (0.26-0.53).
(5) Regla LOC+T (la primaria de C10-C11): solo D2s controla (FDR agrupada 0.00-0.10 a q = 0.05-0.2; la celda
    unilateral en DENSE llega a 0.26 [0.14, 0.41] a q = 0.2) y a costa de casi todo el recall. D1p falla: 0.08-0.19 a
    q = 0.05, 0.18-0.29 a q = 0.1, 0.28-0.37 a q = 0.2 (lote 120), y peor en la celda unilateral (0.34-0.43 en DENSE).
    Motivo: un evento real localizado a <= 5 km pero con el origen mal (> 0.5 s) cuenta como falso, y ningun senuelo
    de moveout lo imita; son 0.2-1.2 atomos por ventana con eventos, segun estrato.
(6) Coste en recall del umbral por senuelos frente al umbral oraculo con la MISMA tasa realizada de falsos: nulo
    con lotes de 120 ventanas (D1p, LOC: perdida entre -0.02 y +0.04 en las lineas agrupadas, intervalos que
    contienen 0); con lotes de 20 ventanas en STD, +0.04 a +0.08 a q <= 0.1 porque el 88% de los lotes no acepta nada.
    El coste real esta en el nivel: a q = 0.1 (LOC, lote 120) el recall es 0.195 en STD y 0.152 en DENSE, frente a
    0.436 y 0.253 en el punto del benchmark (<= 0.25 falsos por ventana), cuya FDP es 0.41 y 0.20. Los puntos de
    operacion del benchmark NO son puntos "sin falsos positivos": un 13-41% (LOC) o un 31-52% (LOC+T) de lo que
    aceptan es falso.
(7) Punto de operacion desplegable equivalente al del benchmark (fpw con D1p, LOC, lote 120): tasa realizada de
    falsos por ventana 0.09 [0.07, 0.12] (STD) y 0.12 [0.09, 0.16] (DENSE) para c = 0.1; 0.21 [0.17, 0.25] y 0.23
    [0.19, 0.28] para c = 0.25; recall 0.355 / 0.181 y 0.453 / 0.245, es decir el del oraculo (0.331 / 0.161 y 0.436 /
    0.253) dentro de +-0.03. Bajo LOC+T la misma regla deja pasar 1.3-3 veces la tasa nominal; D2s la cumple (0.05-
    0.06 para 0.1; 0.12-0.14 para 0.25) con perdida de recall de 0.02-0.05 frente al oraculo.
(8) Deflacion: si. En ventanas con UN evento ya acertado y deflactado, los atomos falsos posteriores a mas de 10 km
    son 1.9 [1.6, 2.4] veces (t80) y 5.0 [3.1, 8.8] veces (t95) los del gemelo sin evento. D3 no lo predice; D1p y D3c
    predicen de mas a t80 (cociente falsos / senuelos 0.67-1.06) porque ven el evento sin deflactar.
(9) bh (p-valor conformal + Benjamini-Hochberg) = tdc, numero por numero. comp (competencia pareada por rango) es
    invalida con cualquier senuelo (D1p, LOC, q = 0.1: 0.23-0.32 frente a 0.07-0.10 de tdc): el atomo r del senuelo
    no es la hipotesis nula del atomo r del blanco.
Veredicto: PROMETEDOR CON CONDICIONES. Hay una calibracion desplegable y barata (una pasada extra del detector con las
estaciones permutadas) que reproduce el punto de operacion del benchmark y controla la FDR a q = 0.2 bajo la regla
LOC en las tres celdas, en ventanas estandar y densas, con clutter fino y fuerte. NO hay garantia teorica (la
intercambiabilidad solo vale sin eventos), NO vale para la regla LOC+T, queda corta en la geometria unilateral a
q <= 0.1, y el objetivo "sin falsos positivos ni falsos negativos" queda lejos: con FDR 0.05 el recall es 0.09-0.10.
Secciones 9-11.
"""

P_BENCH = r"""
----------------------------------------------------------------------------------------------------------------------
3. MUESTRA Y PUNTOS DE OPERACION DEL BENCHMARK (umbral elegido con las etiquetas)
----------------------------------------------------------------------------------------------------------------------
"""
P_BENCH2 = r"""
Lectura: los puntos del benchmark fijan los eventos falsos por ventana, no la proporcion de falsos. En ventanas STD
(0.82 eventos por ventana) 0.25 falsos por ventana equivalen a un 41% (LOC) o 52% (LOC+T) de detecciones falsas; en
DENSE (4 eventos por ventana) a un 20% o 37%. Controlar la FDR a 0.05-0.2 es por tanto mas exigente que el benchmark
en STD y comparable en DENSE bajo LOC.
"""

P_EXCH = r"""
----------------------------------------------------------------------------------------------------------------------
4. INTERCAMBIABILIDAD EN VENTANAS SIN EVENTOS (c13_analysis.txt, primera seccion)
----------------------------------------------------------------------------------------------------------------------
Ventanas de solo clutter: las 8 CLUT de cada red mas las STD sin evento (177-200 ventanas por estrato).
"""
P_EXCH2 = r"""
D1p, D2s, D2r y D3c: ninguna diferencia detectable con el blanco (KS del atomo maximo 0.04-0.10; un solo p por debajo
de 0.05 entre 48, D2r en strong ladder_6 con todos los atomos, p 0.039; cocientes de cola con intervalos que contienen
1). D3: compatible en clutter fino, rechazado en clutter fuerte (p 0.001 en los tres estratos con todos los atomos,
cociente de cola a 0.95 de 1.50-1.75). D1o: rechazado en dos estratos fuertes (p 0.001-0.003) y en thin ladder_4
(p 0.009). D3c, que solo cambia el rango de
los retardos, pasa (p 0.35-0.82; cociente 0.90-1.17): la causa era el borde sin clutter de la traza, no la idea.
Ventanas de solo ruido (24 por estrato): los seis senuelos son compatibles con el blanco (KS 0.08-0.33; una sola de
72 pruebas con p < 0.05, D3c en thin ladder_6, p 0.039; cocientes con intervalos anchos que contienen 1).
Con 24 ventanas la prueba solo detectaria desvios grandes.
"""

P_EV = r"""
----------------------------------------------------------------------------------------------------------------------
5. VENTANAS CON EVENTOS: FALSOS DEL BLANCO FRENTE A SENUELOS Y FRENTE AL GEMELO SIN EVENTOS
----------------------------------------------------------------------------------------------------------------------
"""
P_EV2 = r"""
Lectura de TC. (i) Las detecciones falsas se reparten en tres clases de peso parecido: evento real con el origen mal
(tmiss), lobulo a 5-10 km (near) y lejos (far). Solo la primera desaparece al pasar de LOC+T a LOC. (ii) D3 y D1o
cuentan menos senuelos que falsos (cociente 1.0-2.4 a t80 y mucho mayor en la cola): ademas de perder clutter, el
retardo sobre toda la traza saca de la zona leida gran parte de las llegadas de los eventos. (iii) D1p, D2s, D2r y
D3c cuentan MAS senuelos que falsos a t80 (0.67-1.06): son conservadores en numero, y lo son mas bajo LOC (falsos LOC
/ D1p 0.5-0.9). La distribucion de puntuaciones de los falsos se parece a la de D1p y D3c (KS 0.04-0.19) y no a la de
D3 o D1o (0.20-0.51). (iv) Que el recuento sea conservador a t80 no garantiza la cola: alli decide la FDR (seccion 6).
"""
P_TW = r"""
Gemelos (c13_twin.py): W0 = la misma ventana sin la parte de los eventos (mismo clutter y mismo ruido, muestra a
muestra; normalizacion de W). Atomos por ventana con eventos por encima de t80 | t95:
"""
P_TW2 = r"""
Lectura. Sin eventos, el clutter y el ruido de esas mismas ventanas producen 0.67 (STD) y 0.42 (DENSE) atomos por
ventana por encima de t80 y 0.12 / 0.07 por encima de t95. Con eventos hay 2.3 / 3.9 y 1.4 / 2.9 atomos FALSOS: casi
todo lo falso por encima del nivel del clutter lo crean los eventos (cocientes 3.5 y 9.3 a t80; 12 y 45 a t95), tambien
lejos de ellos (2.3 y 5.2; 6.7 y 24). Pregunta de la deflacion: en ventanas STD con un solo evento que ya fue acertado
y deflactado, los falsos posteriores a mas de 10 km son 1.93 [1.55, 2.39] veces (t80) y 4.97 [3.12, 8.79] veces (t95)
los atomos lejanos del gemelo (thin 2.44 y 6.86; strong 1.41 y 2.78); la comparacion es conservadora porque el gemelo
cuenta todos los rangos. La deflacion (un tensor de momento por atomo, ondicula de autocorrelacion, sin coda) deja
residuo suficiente para sembrar falsos lejos del evento. Un senuelo que no ve la energia de los eventos (D3, D1o) no lo
puede predecir; D1p y D3c la ven antes de deflactar y por eso sobrestiman a umbrales moderados.
"""

P_FDR = r"""
----------------------------------------------------------------------------------------------------------------------
6. CALIBRACION tdc: FDR REALIZADA FRENTE AL NIVEL NOMINAL
----------------------------------------------------------------------------------------------------------------------
Tablas completas (todas las lineas, D1o, q = 0.3 y 0.5, comp, bh, fpw): out/c13_report_tables.txt y
out/c13_analysis.txt. D1o se omite aqui: se comporta como D3 (p. ej. STD agrupado, LOC, lote 120, q = 0.1: 0.42 frente a
0.41). Lineas "both" = los dos regimenes; "ALL" = las tres celdas.
6.1 Regla LOC
"""
P_FDR2 = r"""
6.2 Regla LOC+T
"""
P_FDR3 = r"""
Lectura. Bajo LOC, D1p cumple el nivel en todas las lineas a q = 0.2 con ambos lotes y con las dos formas de intervalo;
a q = 0.05 y 0.1 lo cumple en ladder_4 y ladder_6 y queda en 1.2-2 veces el nominal en unilateral_c40 con lotes de 120
(STD 0.05 y 0.12, sin marca; DENSE 0.10 y 0.15, marcados), no con lotes de 20 (0.00-0.04 y 0.05-0.12). Los lotes de 20
son mas conservadores porque el "+1" del estimador pesa mas. D2r cumple bajo LOC en todas las lineas (0.01-0.05 /
0.04-0.08 / 0.10-0.13 agrupado) salvo DENSE unilateral a q = 0.05 (0.09). D3c se queda corto de forma sistematica
(1.3-2 veces el nominal; a q = 0.2 en DENSE 0.21 [0.18, 0.25], casi nominal). D3 no controla nada. Bajo LOC+T solo
D2s cumple (con la reserva de la celda unilateral en DENSE a q = 0.2: 0.26 [0.14, 0.41]); D1p, D2r y D3c fallan en
DENSE con todos los niveles y la celda unilateral es la peor (D1p 0.34-0.43).
Los intervalos percentiles del bootstrap quedan desplazados hacia arriba (hay estimaciones puntuales por debajo del
extremo inferior de su intervalo, p. ej. 0.02 [0.03, 0.08]); TE2 da la alternativa con errores estandar entre lotes
independientes y lleva a las mismas conclusiones.
6.3 Variantes
"""
P_FDR4 = r"""
bh coincide con tdc en todas las casillas. comp es anticonservadora con todos los senuelos, tambien con los que
controlan con tdc (D1p bajo LOC: 0.12-0.22 a q = 0.05). No se recomienda.
"""

P_REC = r"""
----------------------------------------------------------------------------------------------------------------------
7. RECALL Y COSTE FRENTE AL PUNTO DE OPERACION ORACULO
----------------------------------------------------------------------------------------------------------------------
Columnas de TF: recall | eventos falsos por ventana realizados | recall del umbral oraculo con la misma tasa
realizada | perdida = oraculo - tdc [95%] | fraccion de lotes que aceptan algo.
"""
P_REC2 = r"""
Lectura de TF. Con lotes de 120 la perdida frente al oraculo con la misma tasa de falsos es indistinguible de cero
(LOC, agrupado: STD +0.007 [-0.033, +0.043] a q = 0.1 y +0.003 [-0.039, +0.031] a q = 0.2; DENSE +0.006 [-0.011,
+0.017] y -0.000 [-0.014, +0.017]): el umbral por lote no pierde recall respecto de un umbral global. Con lotes de 20
en STD la perdida es +0.084 [+0.046, +0.146] a q = 0.1 y +0.036 [-0.019, +0.058] a q = 0.2: un lote de 20 ventanas con
unos 16 eventos rara vez junta evidencia para aceptar algo a q <= 0.1 (12% de los lotes; 0% a q = 0.05). En DENSE el
lote de 20 basta (perdida +0.011 y +0.005). El recall absoluto con control de FDR es bajo: LOC, lote 120, q = 0.05 /
0.1 / 0.2 -> STD 0.094 / 0.195 / 0.281, DENSE 0.096 / 0.152 / 0.225; ladder_4 no acepta casi nada a q <= 0.2 (recall
0.00-0.06): con 4 estaciones los blancos verdaderos no se separan de los senuelos.
fpw: tasa realizada de falsos por ventana [95%] (* = por encima de c) | recall | recall del oraculo del benchmark a
<= c | perdida = oraculo - fpw [95%]
"""
P_REC3 = r"""
Lectura de TG. Bajo LOC la estimacion por D1p de los falsos por ventana esta calibrada: ningun estrato queda marcado con
lotes de 120 y las lineas agrupadas dan 0.09-0.12 para c = 0.1 y 0.21-0.23 para c = 0.25, con el recall del oraculo
(perdidas -0.024 [-0.054, +0.007] y -0.017 en STD; -0.020 [-0.048, +0.002] y +0.008 en DENSE; negativo = algo mas de
recall porque la tasa realizada no es identica). La celda unilateral en DENSE es la excepcion (0.18 [0.12, 0.25] para
0.1). Con lotes de 20 la tasa sube a 0.13-0.15 y 0.24-0.26: floor(c x 20) senuelos es una regla gruesa. Bajo LOC+T la
misma regla da 0.17 / 0.32 en STD y 0.32 / 0.51 en DENSE (0.51 / 0.70 en la celda unilateral). D2s bajo LOC+T cumple
con margen (la mitad de la tasa nominal) y pierde 0.048 / 0.050 de recall en STD y 0.016 / 0.034 en DENSE frente al
oraculo; bajo LOC pierde 0.11-0.16. D3c deja pasar 1.0-1.8 veces la tasa nominal bajo LOC y 1.5-4 veces bajo LOC+T.
"""

P_COST = r"""
----------------------------------------------------------------------------------------------------------------------
8. COSTE DE COMPUTO
----------------------------------------------------------------------------------------------------------------------
Segundos por ventana (un proceso; medias por celda, STD | DENSE): pasada del blanco ladder_4 0.08 | 0.12, ladder_6
0.12 | 0.17, unilateral_c40 0.16 | 0.22; las cinco pasadas senuelo juntas 0.41-0.42 | 0.58-0.59, 0.60-0.63 | 0.83-0.87,
0.88-0.89 | 1.17-1.18, es decir cada senuelo cuesta lo que una pasada del detector (el sorteo y las tablas del senuelo
son despreciables). La configuracion propuesta (un senuelo) duplica el coste de deteccion; la calibracion del umbral
es una ordenacion de puntuaciones. Corrida completa de C13: 54 min (blanco + 5 senuelos, 3600 ventanas) + 10 min (D3c)
+ 9 min (gemelos) + 6 min de bootstrap.
"""

P_VERD = r"""
----------------------------------------------------------------------------------------------------------------------
9. JUICIO
----------------------------------------------------------------------------------------------------------------------
Respuestas a las preguntas del encargo.
  (a) Algun senuelo da control valido en TODAS las celdas, incluidas las densas con clutter fuerte? Con la regla LOC,
      si a q = 0.2: D1p (FDR 0.02-0.24 por estrato, 0.14-0.16 agrupada, ninguna marca) y D2r; a q = 0.05-0.1, D1p
      cumple en las celdas ladder y se queda corto en la unilateral con lotes de 120 (hasta 2 veces el nominal). Con la
      regla LOC+T, no: solo D2s, que es valido porque es muy conservador (sus senuelos se alimentan de los eventos
      reales con el tiempo correcto) y no deja recall con tdc (0.00-0.02 a q = 0.1). Dense + strong no es el caso peor;
      lo es la geometria unilateral.
  (b) Que recall se pierde frente al umbral oraculo? Con la misma tasa realizada de falsos, nada medible con lotes de
      120 ventanas y 0.04-0.08 en STD con lotes de 20 a q <= 0.1. Frente al punto del benchmark, el control de FDR a
      0.1 cuesta mas de la mitad del recall en STD (0.195 frente a 0.436) porque el benchmark opera con un 41% de
      falsos; fpw con D1p reproduce el punto del benchmark sin etiquetas y sin perdida bajo LOC.
  (c) La deflacion tras eventos fuertes infla los falsos mas de lo que predicen los senuelos? Infla (1.9 a 5 veces los
      del gemelo, lejos del evento ya deflactado), y un senuelo incoherente sin energia de eventos (D3, D1o) lo
      subestima por un factor grande. D1p y D3c, que ven la energia sin deflactar, no lo subestiman a umbrales
      moderados; en la cola D3c si (FDR 1.3-2 veces el nominal) y D1p solo en la celda unilateral.
Limitaciones que deben viajar con cualquier numero.
  (i)   Ninguna garantia. La intercambiabilidad se demuestra sin eventos; con eventos D1p funciona porque las
        coincidencias parciales de un moveout permutado se parecen a los lobulos reales, y eso es empirico y depende de
        la geometria (falla antes en la red unilateral de 8 estaciones, donde un desarreglo conserva menos estructura).
  (ii)  Depende de la regla de acierto. Bajo LOC+T, la regla primaria de las tareas registradas recientes, el
        estimador deja pasar 1.3-3 veces el nominal. Un evento bien situado con el origen mal no es un "falso
        descubrimiento" en el sentido de una nula sin evento; ningun senuelo de moveout lo va a contar.
  (iii) Una semilla, 12 redes por celda y regimen, 20 ventanas por red. Los lotes de 120 mezclan 6 redes de la misma
        celda; no se probo un lote de 120 ventanas de una sola red.
  (iv)  El detector es MT-FW duro con terminos +SC; no se probo con T3 ni con los terminos JHD de C11, ni con BOMP-PK.
  (v)   Tasa de eventos: STD tiene 0.82 eventos por ventana y DENSE 4. Con sismicidad mas rala el lote necesario para
        aceptar algo crece (con 20 ventanas STD ya no se acepta nada a q = 0.05).
  (vi)  D3 y D1o, tal como se fijaron antes de ver datos, tenian un defecto de construccion; D3c es post hoc y D1o no se
        corrigio. La conclusion "los senuelos incoherentes subestiman" se apoya en D3c (1.3-2 veces el nominal), no en
        los numeros mucho peores de D3 y D1o.
  (vii) Los intervalos percentiles de la FDR estan sesgados hacia arriba; se dan tambien los de error estandar.
  (viii) Multiplicidad sin corregir; la eleccion de D1p como candidato se hace sobre estos mismos datos.

----------------------------------------------------------------------------------------------------------------------
10. NOVEDAD
----------------------------------------------------------------------------------------------------------------------
No es nuevo: la competencia blanco-senuelo (proteomica), los knockoffs, los p-valores conformales con BH, ni el uso de
permutaciones o desplazamientos temporales entre estaciones para estimar coincidencias fortuitas (practica habitual en
asociacion de fases y en detectores por coincidencia). Lo que seria nuevo si se confirma: (a) un punto de operacion
sin etiquetas para un detector de subespacio de tensor de momento sin picado, calibrado con hipotesis de moveout
permutado, con su validez medida celda por celda junto al punto oraculo; (b) la demostracion con gemelos de que los
falsos en ventanas con eventos los generan los eventos (incluido el residuo de la deflacion) y no el clutter, y de
que por eso los sustitutos que destruyen la coherencia subestiman; (c) el limite explicito de la idea bajo una regla
de acierto con tiempo de origen. No es un detector nuevo; es una regla de calibracion y su dominio de validez.

----------------------------------------------------------------------------------------------------------------------
11. CONFIGURACION CONGELADA PROPUESTA Y CRITERIOS PARA UNA CONFIRMACION PRERREGISTRADA (nada de esto se ha corrido)
----------------------------------------------------------------------------------------------------------------------
Solo se propone lo que aqui funciona: regla LOC, senuelo D1p. Para LOC+T no se propone ninguna afirmacion de validez.
Configuracion C13-F (sin constantes libres; codigo c13_lib.py tal como esta):
  detector   MT-FW tope 6, +SC, modo hard (el registrado), kmax como en el benchmark;
  senuelo    D1p, un sorteo por ventana de un flujo propio, mismo kmax, registro = minimo acumulado;
  lote       120 ventanas de la MISMA red y tipo de ventana (en C13 fueron 6 redes x 20);
  reglas     fpw con c = 0.25 (primaria) y c = 0.1; tdc con q = 0.2 (primaria) y q = 0.1;
  acierto    LOC para las afirmaciones; LOC+T se informa sin regla, junto con fpw-D2s como alternativa conservadora.
Datos: dos semillas nunca usadas (a elegir y comprobar con grep antes de registrar; las de los borradores C9-C11 ya
estan comprometidas), las tres celdas, regimenes thin y strong, >= 20 redes por celda y regimen con 120 ventanas STD y
120 DENSE por red de TAGs nuevos; bootstrap jerarquico (redes, luego ventanas; 2000 replicas) y, para la FDR, tambien
media +- error estandar entre redes.
Criterios (cada uno en las DOS semillas, agrupando celdas y regimenes, por tipo de ventana):
  2M-a validez de la tasa: falsos por ventana realizados de fpw-D1p con c = 0.25, cota superior unilateral del 95%
       <= 0.3125 (1.25 c). Piloto: STD 0.207 [0.174, 0.250], DENSE 0.227 [0.192, 0.276].
  2M-b validez de la FDR: tdc-D1p con q = 0.2, cota superior unilateral del 95% <= 0.25. Piloto: STD 0.14 [0.11, 0.23],
       DENSE 0.16 [0.13, 0.20] (percentil); 0.13 [0.09, 0.17] y 0.16 [0.13, 0.19] por error estandar con lotes de 20.
  2M-c coste: recall de fpw-D1p (c = 0.25) menos recall del oraculo a <= 0.25 falsos por ventana, cota inferior
       unilateral del 95% > -0.03. Piloto: +0.017 (STD) y -0.008 (DENSE), con intervalos de +-0.03.
  2M-d guarda por celda (descriptiva con umbral): ninguna celda con tasa realizada puntual > 1.5 c ni FDR puntual
       > 1.5 q; si la celda unilateral lo incumple, la afirmacion se limita a las celdas ladder. Piloto: unilateral
       DENSE 0.28 para c = 0.25 (1.1 c) y 0.18 para q = 0.2; a c = 0.1 y q = 0.1, 0.18 y 0.15 (1.8 c y 1.5 q): por eso los
       niveles primarios son 0.25 y 0.2 y los de 0.1 quedan como secundarios sin regla.
  Afirmacion permitida si a, b, c se cumplen: "un pase senuelo con estaciones permutadas fija sin etiquetas el punto de
  operacion del detector con la tasa de falsos nominal y sin perdida de recall, bajo la regla de localizacion". No se
  registra ninguna afirmacion de "garantia" ni de validez bajo LOC+T.
Potencia aproximada (DE de bootstrap del piloto, 12 redes por estrato, escalada con raiz de redes): 2M-a, DE 0.019-0.021
  en el piloto y 0.015-0.016 con 20 redes: se cumple con probabilidad > 0.95 si la tasa verdadera es <= 0.28 y falla
  con probabilidad > 0.8 si es >= 0.33. 2M-b, DE 0.02-0.03: se cumple con > 0.9 si la FDR verdadera es <= 0.20. 2M-c,
  DE 0.015-0.02: se cumple con > 0.9 si la perdida verdadera es <= 0.
Coste estimado: 2 pasadas por ventana, 0.2-0.45 s por ventana; 2 semillas x 120 redes x 240 ventanas = 57600
ventanas, 4-6 h de CPU.
Valoracion del borrador: registrable solo si la tesis acepta LOC como regla de la afirmacion. Si la regla primaria
sigue siendo LOC+T, este resultado sirve como diagnostico (seccion 5) y como alternativa conservadora (fpw-D2s), no
como contribucion.

----------------------------------------------------------------------------------------------------------------------
12. ARCHIVOS (todos en /home/claude/work/explore_contrib; salidas en out/)
----------------------------------------------------------------------------------------------------------------------
  c13_lib.py            senuelos D1p, D1o, D2s, D2r, D3; reglas tdc, bh, comp; constantes y TAGs
  c13_run.py            corrida principal -> out/c13_parts/REGIMEN_CELDA_J.json (72 partes; log out/c13_run.log)
  c13_run_d3c.py        senuelo corregido D3c -> out/c13_d3c_parts/ (log out/c13_d3c.log)
  c13_twin.py           gemelos sin eventos -> out/c13_twin_parts/ (log out/c13_twin.log);
                        out/c13_twin_parts_BUGGY_v1 = primera version invalida, conservada, no usada
  c13_analyze.py        intercambiabilidad, ventanas con eventos, gemelos, calibracion y bootstrap ->
                        out/c13_analysis.txt, out/c13_summary.json (regla fpw incluida)
  c13_tables.py         tablas estrechas -> out/c13_report_tables.txt
  c13_report.py, c13_report_text.py   montaje de este informe
  c13_queue.sh          cola posterior a la corrida (un proceso a la vez)
  out/c13_analysis_v1_5decoys.txt, out/c13_report_tables_v1_5decoys.txt, out/c13_summary_v1_5decoys.json: analisis
                        previo a D3c y con los gemelos invalidos (solo historico); out/c13_analysis_interim4.txt: mirada
                        intermedia con 4 redes
  No se modifico ningun archivo registrado, ninguna parte 2A-2L, ningun archivo anterior de explore_contrib ni nada en
  /mnt/user-data/outputs.
"""

BOTH = lambda x: " both " in x          # noqa: E731

PARTS = [
    P_HEAD, P_SUM, P_DESIGN, P_DECL,
    P_BENCH, ("section", "TA"), P_BENCH2,
    P_EXCH, ("sub", "TB", "  clutter-only"), P_EXCH2,
    P_EV, ("section", "TC"), P_EV2, P_TW, ("td", "STD both ALL"), ("td", "DENSE both ALL"), P_TW2,
    P_FDR,
    ("sub", "TE", "  rule LOC, batch 120 windows, q = 0.05", BOTH),
    ("sub", "TE", "  rule LOC, batch 120 windows, q = 0.1"),
    ("sub", "TE", "  rule LOC, batch 120 windows, q = 0.2"),
    ("sub", "TE", "  rule LOC, batch 20 windows, q = 0.1", BOTH),
    "TE2 (intervalos por error estandar entre lotes independientes de 20 ventanas):",
    ("sub", "TE2", "  rule LOC, batch 20 windows, q = 0.1", BOTH),
    ("sub", "TE2", "  rule LOC, batch 20 windows, q = 0.2", BOTH),
    P_FDR2,
    ("sub", "TE", "  rule LOC+T, batch 120 windows, q = 0.05", BOTH),
    ("sub", "TE", "  rule LOC+T, batch 120 windows, q = 0.1"),
    ("sub", "TE", "  rule LOC+T, batch 120 windows, q = 0.2", BOTH),
    ("sub", "TE2", "  rule LOC+T, batch 20 windows, q = 0.1", BOTH),
    P_FDR3,
    ("sub", "TH", "  STD both ALL, rule LOC, batch 120"),
    ("sub", "TH", "  DENSE both ALL, rule LOC, batch 120"),
    ("sub", "TH", "  DENSE both ALL, rule LOC+T, batch 120"),
    P_FDR4,
    P_REC,
    ("sub", "TF", "  decoy D1p, rule LOC, batch 120"),
    ("sub", "TF", "  decoy D1p, rule LOC, batch 20", BOTH),
    ("sub", "TF", "  decoy D1p, rule LOC+T, batch 120", BOTH),
    "TF0 (D1p, q = 0.05 y q = 0.5):",
    ("sub", "TF0", "  rule LOC, batch 120", BOTH),
    P_REC2,
    ("sub", "TG", "  decoy D1p, rule LOC, batch 120"),
    ("sub", "TG", "  decoy D1p, rule LOC, batch 20", BOTH),
    ("sub", "TG", "  decoy D1p, rule LOC+T, batch 120", BOTH),
    ("sub", "TG", "  decoy D2s, rule LOC+T, batch 120", BOTH),
    ("sub", "TG", "  decoy D3c, rule LOC, batch 120", BOTH),
    P_REC3, P_COST, P_VERD,
]
