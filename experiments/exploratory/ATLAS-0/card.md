# Experiment card: ATLAS-0. The failure atlas of published detection pipelines

Stage: S0 (draft of 18 October 2026). It must be completed and anchored (ADR-0004) before any comparator of the atlas
is run, because the gate of phase 4 counts wins over "the regimes of the atlas", and regimes chosen after seeing the
atlas would make that gate circular.
Owner decisions: pending, see "Open choices".

## Question
In which measurable regimes does each published detection pipeline lose catalogued events, and how many?

## Mechanism tested
The failure of each pipeline stage: picker (missed or merged picks against arrival separation and signal-to-noise
ratio), associator (merged and split events, minimum-pick thresholds), velocity model.

## Hypothesis and refutation
H: there is at least one regime of the grid below in which every core pipeline loses at least 20 per cent of the
events that the oracle injection says are detectable. Refuted if no cell meets it with its 95 per cent interval.

## Regime grid (fixed here, not after the atlas)
Each injected or catalogued event is assigned to one cell of:
- separation to the nearest other event's arrival at the median station: < 0.5 s, 0.5 to 1 s, 1 to 3 s, > 3 s;
- network signal-to-noise ratio of the event (median over stations, P window, dB): < 0, 0 to 6, 6 to 12, > 12;
- number of stations within 50 km: <= 6, 7 to 15, > 15.
That gives 48 cells. Cells with fewer than 100 events (P6-power) are reported but not counted.

## Counting rule for "majority of regimes" (used by the gate of phase 4)
A method wins a cell when its recall under LOC+T (ADR-0007) at the cell's budget exceeds the best other method's
recall with a paired 95 per cent interval above zero. "Majority" is a weighted majority of counted cells, each cell
weighted by its number of events. Ties and cells with intervals covering zero count as not won.

## False-detection budget
Primary: 30 false detections per hour of record, as in R1; secondary: 10 and 100. Every method also reports its
deployable operating point, fixed on development.

## Truth design
1. Semi-synthetic: recorded events (waveforms of catalogued events at each station) injected into real continuous
   noise at controlled separations and scales. The event library used for injection is disjoint, by event and by time
   block, from any library a method may use for templates or deflation.
2. Synthetic generator v2, independent of the detector's forward model: another wavelet family, a layered or 3-D
   medium and a radiation operator that is not the detector's far-field operator.
3. Real catalogues as secondary evidence only. A catalogue produced by template matching is not used to evaluate
   template matching.

## Methods and comparators
Core set: PhaseNet and PhaseNO with GaMMA and PyOcto; EQNet; template matching (EQcorrscan); BOMP-PK and MT-FW.
Extended set if time allows: EQTransformer, REAL, QuakeMigrate, MALMI. PhaseNet peak distance and threshold, GaMMA and
PyOcto minimum picks and every associator's velocity model are tuned on development with the same budget for every
method (protocol, Section 3).

## Data
To be listed by data card before the run. At least two networks; development partitions only (phase 3 is not
confirmatory). Data cards of the held datasets: `data/cards/`.

## Gate
The harness reproduces published numbers of at least two comparators on their own example data, within 5 per cent
of their recall, before any atlas cell is read.

## Open choices
- The grid edges above are a default; changing them is allowed only before the card is anchored.
- Whether the minimum of 100 events per counted cell can be met depends on the data (ADR-0005).

## Cost and risks
PhaseNO and EQNet inference on CPU is slow for long records; training is not part of the atlas. The main risk is too
few events in the overlapping cells; the semi-synthetic injection is what fills them.
