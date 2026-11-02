# Prerregistro de la fase 6 en OSF: borrador campo por campo

BORRADOR, no registrado. 2 de noviembre de 2026. Las instrucciones están en español; el texto para pegar en OSF está en
inglés, como el resto de los registros del repositorio (`preregistration/drafts/`). Este archivo no cambia ninguna
regla: copia en los campos de OSF lo que ya fijan `docs/decisions/ADR-0008-single-beat-rule.md`,
`docs/decisions/ADR-0004-external-anchor.md`, `docs/protocol/PHASE6_SAP.md` (v0.1) y los borradores
`preregistration/drafts/taskR2_registration_DRAFT.md` y `taskR3_registration_DRAFT.md`. Si alguno de ellos cambia,
este borrador se actualiza a partir de él, no al revés.

## 0. Lo que hay que saber antes de abrir OSF

1. **Plantilla recomendada: "OSF Preregistration"** (la plantilla general de OSF, sucesora del Preregistration
   Challenge). Es la que tiene campos separados para hipótesis, datos existentes, tamaño de muestra, regla de parada,
   variables, criterios de inferencia, exclusiones y datos faltantes, que son justo las piezas del SAP. Alternativas
   descartadas: "Secondary Data Preregistration" (pensada para bases de datos de terceros ya analizadas; R3 crea sus
   datos y R2 no ha leído las formas de onda), "Open-Ended Registration" (un solo campo libre: sirve solo para el
   registro opcional del punto 4), "Preregistration in Social Psychology" y las de ensayos clínicos (no aplican).
2. **Dos registros, no uno.** R2 (datos reales, red diezmada) y R3 (semisintético, redes no vistas) responden
   preguntas distintas, cada una a 0.025 (SAP §1). Cada uno es su propio registro OSF dentro del mismo proyecto.
3. **Hoy se puede crear el proyecto, pero todavía no se pueden enviar R2 ni R3.** Un registro OSF es inmutable y no se
   puede borrar; enviar uno con huecos lo deja inútil. Faltan, en este orden (plan maestro, fase 6):
   - G4 cumplido en desarrollo con la regla de ADR-0008 (fase 4); de ahí salen la predicción y los rangos de
     `sample_size.py estimate`;
   - la tarjeta de datos sellada de las secuencias reservadas (T3.R4) dimensionada por T6.2/T6.3;
   - sigma_t de la verdad medido en T3.V (decide si R2 puede probar T);
   - el nivel a2 de H2 del chequeo de diseño (`sample_size.py calibrate`) y el número de bloques (`solve`);
   - el código congelado con su tag `reg/R2/freeze-1` o `reg/R3/freeze-1`.
   Cada hueco está marcado abajo como ⟨…⟩ con la fuente que lo llena. Ninguno se puede llenar con datos reservados.
4. **Opcional, se puede hacer ya: registrar el plan de análisis.** Un registro "Open-Ended Registration" con
   ADR-0008 y el SAP v0.1 adjuntos fija con fecha y DOI de un tercero que la regla de victoria existía antes de los
   resultados de desarrollo de G4. No obliga a ningún número y, si el SAP cambia, el cambio queda como enmienda
   visible. Mi recomendación es hacerlo, pero es decisión del autor. Texto propuesto en la sección 5.
5. **Orden de anclaje (ADR-0004).** Para R2 y R3: commit del registro lleno → tag anotado y Release
   `reg/R2/registration-1` → `ots stamp` del SHA-256 → registro OSF con el SHA-256, el commit y la URL de la Release
   en el texto → solo después, la corrida. El DOI de OSF se anota en el log de resultados.

## 1. Pasos en OSF (los hace el autor)

1. Entrar en https://osf.io con la cuenta ligada a ORCID 0009-0005-4252-1263 (perfil https://osf.io/geauh/).
2. **Crear el proyecto**: "My Projects" → "Create new project". Título:
   `Timing-tolerant moment-tensor detection of overlapping earthquakes`. Región de almacenamiento: la que ofrezca por
   defecto (no importa para el registro). Dejarlo **privado** hasta decidir lo contrario; un registro puede hacerse
   desde un proyecto privado.
3. En el proyecto, llenar los metadatos (sección 2 de este archivo): descripción, licencia CC BY 4.0, sujetos,
   etiquetas, institución. Añadir como contribuidor al director (Uram Aníbal Sosa Aguirre) si tiene cuenta OSF; si
   no, se puede añadir sin cuenta por correo o dejarlo fuera y nombrarlo en la descripción. Todo contribuidor
   administrador debe aprobar cada registro en 48 h o se aprueba solo.
4. **Add-ons → GitHub**: conectar `juanjgc0609/timing-tolerant-mt-detection`. Al registrar, OSF archiva una copia
   del repositorio enlazado tal como esté en ese momento; por eso el registro se envía después de empujar el tag.
5. **Files → OSF Storage**: subir, para cada registro, el archivo de registro lleno y su `.ots`, `PHASE6_SAP.md`,
   `ADR-0008-single-beat-rule.md`, `CR04_PREVALENCE.md`, la tarjeta de datos sellada y la salida de
   `sample_size.py` usada. Así el registro no depende de que GitHub siga en pie.
6. **Registrar**: pestaña "Registrations" → "New registration" → plantilla "OSF Preregistration". Pegar los
   campos de la sección 3 (R2) o 4 (R3). En "Files" de cada pregunta que lo permita, adjuntar los archivos del paso 5.
7. **Embargo**: OSF ofrece registro público inmediato o embargo de hasta 4 años. Recomiendo **público inmediato**
   para R3 (sus semillas se derivan de la Release y no revelan nada). Para R2 también público, salvo que la licencia
   de los datos reservados impida nombrar las secuencias; en ese caso, embargo hasta la publicación del resultado.
   El embargo no cambia la fecha del registro.
8. Revisar la vista previa, enviar, y copiar el DOI (formato `10.17605/OSF.IO/XXXXX`) al log de resultados y a la
   tabla de la sección 6.

## 2. Metadatos del proyecto y de cada registro

| Campo OSF | Texto |
|:--|:--|
| Title (proyecto) | Timing-tolerant moment-tensor detection of overlapping earthquakes |
| Title (R2) | Confirmatory test R2: overlapping-event detection on decimated real networks |
| Title (R3) | Confirmatory test R3: overlapping-event detection with recorded events injected into unseen networks' noise |
| Description (proyecto) | Doctoral research project (Universidad Icesi, director Uram Aníbal Sosa Aguirre) on detecting and locating earthquakes whose arrivals overlap within 0.1 to 1 s, where pick-and-associate methods lose events. Code, protocol, decision records and the sealed registration log are in the linked GitHub repository; confirmatory tests of phase 6 are preregistered here (ADR-0004). |
| Description (R2/R3) | One paragraph: the "Question" of the registration (sección 3/4, campo Hypotheses, primera frase) plus "Analysis rule: ADR-0008; analysis plan: PHASE6_SAP v⟨versión⟩; anchored by tag and Release reg/R⟨2/3⟩/registration-1 (commit ⟨hash⟩, SHA-256 ⟨hash⟩)." |
| Contributors | Juan José Gordillo Córdoba (ORCID 0009-0005-4252-1263); Uram Aníbal Sosa Aguirre (director) si acepta figurar |
| Affiliated institutions | Universidad Icesi (si OSF la lista; si no, en la descripción) |
| License | CC BY 4.0 (texto y documentos). El código del repositorio enlazado es MIT; decirlo en la descripción. |
| Subjects | Physical Sciences and Mathematics → Earth Sciences → Geophysics and Seismology; Physical Sciences and Mathematics → Statistics and Probability |
| Tags | seismology; earthquake detection; overlapping events; moment tensor; non-inferiority; preregistration; benchmark |
| Category | Project (proyecto); Analysis o Methods and Measures (componentes, si se crean) |

## 3. Registro R2 (datos reales), campo por campo

Fuente: `preregistration/drafts/taskR2_registration_DRAFT.md` y `docs/protocol/PHASE6_SAP.md`. Lo que va entre
⟨⟩ se llena antes del anclaje; nada cambia después de leer una forma de onda reservada.

### Study Information

**Hypotheses**

> Question. On real sequences from networks never used in development, does the model detect more catalogued events
> than every tuned comparator in the target regime T (arrival separation to the nearest other event at the median
> station in [0.1, 1) s) and no fewer, by more than 0.05 of recall, in the complement C (every other event), at 30
> false detections per hour?
>
> H1 (tested only if T is testable on this truth, see Sampling Plan): for every comparator k, recall in T of the
> model minus recall in T of k is greater than 0. Directional, one-sided.
> H2 (tested only if H1 is rejected, or as the primary hypothesis if T is not testable): for every comparator k,
> recall in C of the model minus recall in C of k is greater than -0.05 (non-inferiority, margin 0.05).
> Secondary hypotheses (one Holm family at 0.025, tested only if H1 and H2 are rejected): superiority in [0.1, 0.5) s;
> superiority in [0.5, 1) s; non-inferiority in C with margin 0.03; the primary rule at 10 and at 100 false detections
> per hour; the primary rule at the deployable threshold fixed on development.
> Prediction (from development, phase-4 gate G4, report ⟨ruta⟩): ⟨efecto previsto en T y en C⟩.
> Rationale: Appendix AB (identifiability of overlapping sources) predicts a gain in T and none below 0.1 s, which is
> below the resolvable limit (Cor. A.7).

### Design Plan

**Study type**: Observational Study (methods are compared on recorded data; nothing is assigned).

**Blinding**: No blinding is involved in this study (no human subjects).

**Is there any additional blinding in this study?**

> The evaluation sequences were chosen and sealed (data card with file list and SHA-256, anchored ⟨tag⟩) before any of
> their waveforms was read. Only the catalogue and station metadata were read, to compute eligibility, regimes and
> the expected number of targets. All method code, configurations, tuning ledgers and the analysis script are frozen
> (tag reg/R2/freeze-1) before the run. The run is single.

**Study design**

> Paired comparison of detection methods on the same records. Methods: the model ⟨tag, versión⟩ and the core
> comparator set of ATLAS-0 (registration reg/ATLAS-0/registration-1, commit 7bafb7a, amendments 1 to 3, the last at
> commit 83c66e4): PhaseNet and PhaseNO pickers with GaMMA and PyOcto associators, EQNet, template matching where the
> truth is not its own product, BOMP-PK and MT-FW; plus the best comparator re-run with the model's picker if that
> picker is fine-tuned. Each method gets the same tuning budget: 60 configurations by the same random search, seed
> list and development objective, every trial in a tuning ledger (ledger ⟨ruta, SHA-256⟩). All methods receive the
> same waveforms, station list, velocity model and typed station terms. Evaluation on station subsets of 4 to 12
> stations of the dense or nodal network ⟨nombre⟩, excluding the stations used to build the truth (network
> decimation), subsets drawn by ⟨regla⟩ before the run. Replication unit: consecutive 1-hour time blocks within each
> sequence (T/30 when a sequence holds fewer than 30 hours).

**Randomization**: No randomization (observational). Station subsets are drawn by ⟨regla y semilla⟩ fixed before the
run.

### Sampling Plan

**Existing Data**: Registration prior to analysis of the data.

**Explanation of existing data**

> The waveforms exist in public archives but no waveform of the reserved sequences has been read by any person or
> method (data card ⟨tag⟩, SHA-256 of every file). The reviewed catalogue and station metadata have been read, as the
> protocol allows, to compute eligibility, regimes and expected numbers of targets. The data card's anchor time
> precedes this registration; its SHA-256 list lets anyone check that the evaluated files are the sealed ones.

**Data collection procedures**

> Sequences ⟨lista⟩ of networks ⟨lista⟩, periods ⟨UTC⟩, downloaded from ⟨centros de datos⟩ by the frozen script
> ⟨ruta⟩. Development partition: records ending at least 1 h before the evaluation period (causal split): anything a
> method estimates from events of the test network (station terms, template dictionary, moment-tensor prior,
> background statistics) uses only those records. Evaluation partition: ⟨límites UTC exactos, márgenes incluidos⟩.
> Pretrained weights: provenance table ⟨ruta⟩; a model whose training data overlap the test network and period is
> replaced or its overlap is declared.

**Sample size**

> ⟨B⟩ blocks of 1 h with ⟨n_T⟩ targets in T and ⟨n_C⟩ in C expected from the catalogue, ⟨h⟩ hours of record.

**Sample size rationale**

> Simulated with experiments/exploratory/P6-SAP/sample_size.py, which implements this exact rule, for this design:
> recall ⟨p0⟩, method correlation ⟨rho⟩, intra-block correlation ⟨icc⟩ and targets per block ⟨m, k⟩ estimated on
> development hits (sample_size.py estimate, output ⟨ruta⟩), share of targets in T ⟨f⟩, planning effect +0.10 in T
> and 0 in C. Power of the chain H1 then H2: ⟨p⟩; simulated size: ⟨s⟩ (output ⟨ruta⟩). For reference, with f = 0.3,
> rho 0.7 and icc 0.1 the chain needs about 525 blocks (470 targets in T, 1,570 in total) at +0.10; non-inferiority at
> margin 0.05 needs about 770 independent targets in C per comparator.

**Stopping rule**

> Single run on the full reserved data; no data are added after any result is read. If, before the run, the eligible
> targets fall below the number that gives power 0.50, the run does not start and an amendment is required. One
> interim look, after the first reserved sequence, for futility only and non-binding: if the conditional power of H1
> (of H2 when T is not testable) at the registered effect is below 0.10, the author may stop and report the test as
> stopped for futility. If access to data is lost, the test is reported with the blocks run and its achieved power.

### Variables

**Manipulated variables**: None (observational). The compared factor is the detection method.

**Measured variables**

> For each target event and method: hit or miss under the hit rule LOC+T-V (ADR-0008 item 2): one-to-one matching in
> decreasing score order, epicentre within 5 km and origin within 0.5 s, and a detection is eligible for an event only
> if its normalised cost max(distance / 5 km, |origin error| / 0.5 s) to that event is at most half its cost to every
> other catalogued event. Per method: the score of each detection and the number of false detections per hour.
> Per target, computed from the catalogue only, never from method outputs: arrival separation to the nearest event at
> the median station, signal-to-noise ratio, station subset, distance, K_mec/E.

**Indices**

> Primary endpoint: recall of the targets of a regime at the operating point where each method has at most
> floor(30 x hours) false detections over the evaluated record (exact FROC; tied scores enter together), pooled over
> the regime's cells with weights equal to each cell's share in real dense sequences (CR-04 prevalence, table ⟨ruta⟩)
> divided by its share among the evaluated targets. Per block j and comparator k: d_j = weighted sum over the block's
> targets of (hit by model - hit by k); b_j, c_j = targets hit only by the model, only by k; n_j = targets.

### Analysis Plan

**Statistical models**

> H1: block sign-flip test on d_j (exact randomisation form of the clustered McNemar test, Durkalski et al. 2003),
> exact when at most 16 blocks are non-zero, otherwise 100,000 random sign vectors from seed 20261023; the clustered
> McNemar z is reported next to it. H2 and the margin-0.03 secondary: clustered score test on the per-block
> discordant counts, statistic sum_j (b_j - c_j) + margin x n_j, variance under the restricted null (cluster-robust,
> floored at Tango's 1998 multinomial variance at the restricted MLE), Student t reference with Satterthwaite degrees
> of freedom (Bell and McCaffrey 2002). Intervals: stratified block bootstrap, 10,000 replicates, seed 20261023,
> strata = sequence, operating points and best comparator recomputed per replicate, percentile 95 %. Code:
> src/ttmt/stats.py (signflip_pvalue, ni_score_pvalue, block_sums, block_counts, iut, holm, block_bootstrap,
> budget_check), frozen analysis script ⟨ruta, SHA-256⟩.

**Transformations**: None beyond the prevalence weights above.

**Inference criteria**

> One-sided. Beating the comparators means beating each of them: intersection-union test, p-value = largest pairwise
> p-value (Berger 1982), no correction across comparators. Fixed sequence: H1 at 0.025, then H2 at level a2 = ⟨a2⟩
> fixed by the design check (sample_size.py calibrate over the registered range ⟨rangos de p0, m, k, rho, icc⟩: the
> largest of 0.025, 0.02, 0.015, 0.01, 0.005 whose simulated size at the least favourable configuration, plus 1.645
> Monte Carlo SE, is at most 0.025). H2 needs at least 30 blocks. Secondary Holm family at 0.025 after H1 and H2, the
> margin-0.03 threshold scaled by a2 / 0.025. Negative control, mandatory: the duplicator of the best comparator (each
> detection emitted twice at origin +- half the cell separation); if it gains recall in T over the comparator it
> duplicates (one-sided 0.025), the metric is declared broken and no beat claim is made. Deployable-threshold
> secondary: claimed only for a method whose realised false-alarm rate passes the quasi-Poisson score test against
> 1.2 x budget at one-sided 0.025, reported stratified by blocks with and without a parent event.
> Interpretation: H1 and H2 met: superior in T, which holds ⟨Y⟩ % of the events of dense sequences, and not inferior
> by more than 0.05 elsewhere, on this real data. H2 met alone: not inferior on real data; the T claim rests on R3.
> H2 not met: worse than at least one tuned comparator outside T, reported with the discordant targets.

**Data exclusion**

> Decided before the run by the frozen eligibility script ⟨ruta, SHA-256⟩. A station-hour is unusable with more than
> 5 % gaps, flat fill (identical samples for at least 1 s, replaced by gaps) or clipped values, or with missing
> response or orientation. Targets: catalogue events with M >= Mc(t) (maximum curvature + 0.2, raised after large
> shocks by Helmstetter et al. 2006) and at least 4 usable stations of the subset over [origin - 10 s, origin + 60 s].
> Other events are excluded events: a detection matched to them is neutral. A block enters only if at least 50 % of its
> station-hours are usable; excluded blocks are listed with their cause before the run. Truth precision: the truth's
> origin-time standard deviation sigma_t, measured on T3.V, is ⟨valor⟩; a pair enters T only if its separation is at
> least 9 sigma_t (⟨n⟩ pairs qualify); if none qualifies, R2 tests only H2.

**Missing data**

> A method that fails on a block (crash, timeout) scores no detection there (counted as misses); a sensitivity
> analysis removes such blocks for every method. Unusable station-hours are gaps, never zeros.

**Exploratory analysis**

> Descriptive, not tested: every atlas cell with its estimate and interval (cells are never counted), the
> epicentre-only rule, duplicates of isolated events per method, origin and location errors, per-sequence estimates.
> Sensitivity: a detection that matches no catalogued event is neutral when it matches, by LOC+T, an event of the
> independent enhanced catalogue ⟨nombre⟩ (never a template-matching catalogue to score template matching).

### Other

> Analysis rule: docs/decisions/ADR-0008-single-beat-rule.md. Analysis plan: docs/protocol/PHASE6_SAP.md v⟨versión⟩,
> Sections 3 to 10 part of this registration. Registration file ⟨ruta⟩, SHA-256 ⟨hash⟩, commit ⟨hash⟩, tag and GitHub
> Release reg/R2/registration-1 (created_at ⟨UTC⟩), OpenTimestamps proof ⟨ruta .ots⟩. Freeze: reg/R2/freeze-1.
> Deviations known before the result are logged before it is read. Verification: reimplementation by a separate AI
> session (so described) plus the hit or miss of ten discordant targets re-derived by a person, chosen by the seed of
> ADR-0004 item 4. Declared limits: ⟨red, región, periodo, rango de magnitud⟩; truth from one network; decimation is
> not an independent network. Code MIT, text CC BY 4.0.

## 4. Registro R3 (semisintético), solo lo que cambia frente a R2

Fuente: `preregistration/drafts/taskR3_registration_DRAFT.md`. Los campos no listados se copian de la sección 3
cambiando R2 por R3.

**Hypotheses**: misma estructura, con esta pregunta y sin la condición "only if T is testable":

> Question. When recorded events of networks never used in development are injected into those networks' own
> continuous noise at controlled separations, does the model detect more events than every tuned comparator at
> arrival separations of 0.1 to 1 s and no fewer, by more than 0.05 of recall, elsewhere, at 30 false detections per
> hour?

**Study type**: Experiment (the separations, signal-to-noise ratios and station subsets are set by a design grid and
random seeds).

**Additional blinding**: igual que R2, más: "the injection seeds are the first 8 bytes of SHA-256(⟨commit del tag⟩ ||
⟨created_at de la Release⟩) and could not be known before the registration was public (ADR-0004 item 4)."

**Study design**

> Same methods, tuning budget and shared inputs as R2. Reserved networks ⟨lista⟩ (T3.R4). Noise library: ⟨horas⟩ of
> continuous record with no catalogued event above ⟨M⟩ within ⟨km⟩. Event library: ⟨n⟩ recorded events, disjoint by
> event and by time block from every library a method may use for templates, deflation or priors, with no injected
> event correlating above ⟨c⟩ with a library event (CR-06). Design grid: separation in 5 bins, signal-to-noise ratio,
> station subsets of 4, 6, 8 and 12, distance to the calibration centroid, K_mec/E as covariate. One block = one
> network realisation (one noise block, one station subset, its injected events).

**Randomization**: "Injection times, separations within bins, event draws and station subsets are drawn from the
seed above by ⟨script, SHA-256⟩."

**Existing Data**: Registration prior to creation of data.

**Explanation of existing data**: "The injected data set is created only after the anchor, from seeds derived from
it. Its sources (continuous noise and recorded events of the reserved networks) exist in public archives; none of
their waveforms has been read (data card ⟨tag⟩)."

**Sample size / rationale**: "⟨B⟩ blocks from sample_size.py solve with rho ⟨r⟩, icc ⟨i⟩, f ⟨f⟩ and +0.05 in T
(output ⟨ruta⟩); power of the chain ⟨p⟩, size ⟨s⟩." Referencia: con f = 0.3, rho 0.7 e icc 0.1, +0.05 pide unos
1,100 bloques (992 objetivos en T, 3,311 en total).

**Stopping rule**: "Single run on all ⟨B⟩ blocks; no interim look; the run does not start if the realised design
gives power below 0.50 (amendment required)."

**Manipulated variables**: "Arrival separation (5 bins), signal-to-noise ratio, station subset size (4, 6, 8, 12),
distance to the calibration centroid, set by the design grid; the method compared."

**Statistical models**: igual que R2, con bootstrap estratificado por red.

**Inference criteria**: igual que R2, más el chequeo de e-valores: "for the model's e-value rule, the rate of
S_mix >= log(1/p) on hypotheses placed at 3 to 15 km and 1 to 4 s from injected events is reported, and the M7
guarantee is declared refuted on those hypotheses when that rate exceeds 2p." Interpretación: "H1 and H2 met:
superior in T and not inferior elsewhere under known truth on unseen networks. H1 not met: no demonstrated gain where
the theory predicts one; the cells are described. Duplicator control positive: no claim."

**Data exclusion**: sin Mc(t) ni sigma_t (la verdad es la inyección); mismas reglas de estación-hora y bloque.

**Other**: tag `reg/R3/registration-1`, freeze `reg/R3/freeze-1`; "Declared limits: injection adds signals linearly
to real noise; no wave interaction; events of one network are injected into its own noise."

## 5. Registro opcional del plan de análisis (se puede enviar hoy)

Plantilla "Open-Ended Registration". Título: `Analysis rule and statistical analysis plan for the phase-6
confirmatory tests (before development results)`. Campo único:

> This registration fixes, before any phase-4 development result, the rule by which the confirmatory tests R2 and R3
> will decide whether the model beats the comparators: docs/decisions/ADR-0008-single-beat-rule.md (accepted
> 24 October 2026) and docs/protocol/PHASE6_SAP.md v0.1 (23 October 2026), attached, at commit ⟨hash de main⟩ of
> https://github.com/juanjgc0609/timing-tolerant-mt-detection. In short: target regime T = arrival separation
> [0.1, 1) s at the median station; hit rule LOC+T-V; endpoint recall at 30 false detections per hour, weighted to
> real prevalence; H1 superiority in T (block sign-flip), then H2 non-inferiority by 0.05 in the complement (clustered
> score test at the level fixed by a design check), each an intersection-union over every comparator, one-sided
> 0.025; mandatory duplicator negative control; equal tuning budget of 60 configurations per method. The exploratory
> failure atlas that motivates the tests is anchored separately (tag reg/ATLAS-0/registration-1, commit 7bafb7a;
> amendments 1 to 3, the last at 83c66e4). This is not the registration of R2 or R3: each will be registered with its
> data, sample size and level before its run, and any change to the rule will be an amendment that cites this one.

Adjuntar ADR-0008, PHASE6_SAP.md, CR04_PREVALENCE.md y `experiments/exploratory/P6-SAP/sample_size_out.txt`.

## 6. Registro de anclajes (se llena al enviar)

| Registro | Plantilla | Tag y Release | SHA-256 | DOI OSF | Fecha |
|:--|:--|:--|:--|:--|:--|
| Plan de análisis (opcional) | Open-Ended Registration | — (commit de main) | — | | |
| R2 | OSF Preregistration | reg/R2/registration-1 | | | |
| R3 | OSF Preregistration | reg/R3/registration-1 | | | |
