# Detecting and locating earthquakes without travel times: a literature scan

Date: 2 Oct 2026. Scope: 2000-2026, eight method families, with a ranked shortlist for the thesis.
Status: literature scan, not a systematic review. No experiment was run for this document.

## 0. How to read this document

Verification. Every entry in the reference list was opened during this scan on at least one page
(Crossref, OpenAlex or Semantic Scholar record, publisher page, arXiv, or an institutional
repository). Each entry carries a tag:

- [A] abstract or full page read; statements about method, accuracy and limits come from it.
- [M] only bibliographic metadata was read (title, authors, journal, year, DOI). For these, what
  the text says about the paper is limited to what the title supports, plus general background
  that is marked as such.

Volume and page numbers are given only where the opened page showed them.

Sentences that begin with "Background:" are standard seismology or my own reasoning. They were
not checked against a source in this scan and should be treated as claims to verify.

Pages that could not be read (so nothing is cited from them): AGU/Wiley article pages (HTTP 403),
Science.org (403), GeoScienceWorld (403), IOPscience (robots.txt), HAL (access denied), NASA ADS
(robots), GEOMAR OceanRep (JavaScript wall), earth-prints (empty page), and several Crossref,
OpenAlex and Semantic Scholar queries that returned HTTP 429. Papers I intended to include but
could not open, and therefore left out: Lomax, Michelini and Jozinovic 2019 (single-station CNN),
Cristiano et al. 2016 (P polarization at the Graefenberg array), Corciulo et al. 2012 and Wang,
Templeton and Harris 2015 (matched field), Ito, Jin and Zou 2012 (direct sampling), Guilhem and
Dreger 2011 and Tsuruoka et al. 2009 (continuous moment-tensor scanning), Hu, Kian and Zhao 2020
(uniqueness for wave-equation sources). They are worth opening, but they are not evidence here.

The two questions asked of each family:

- Does it truly avoid travel times, or does it hide them (for example in the phase of a Green's
  function)?
- Can it add information independent of arrival times (polarization direction, amplitude decay,
  radiation pattern, wavefield gradients)?

## 1. Summary table

| Family | Observable | Avoids travel times? | Independent information | Works on 4-8 station 3C network? |
|---|---|---|---|---|
| 1 Polarization | P (and S) particle-motion direction | Yes for direction; distance still needs S-P or amplitude | Yes: azimuth, and incidence angle if near-surface velocity is known | Yes, station by station; reliability varies by site |
| 2 Amplitude (ASL, ratios, S/P) | Band-limited RMS or envelope amplitude per station | Yes (coarse window alignment only) | Yes: distance decay, radiation pattern | Yes, but site terms dominate the error budget |
| 3 Time reversal, RTM, adjoint, waveform inversion | Full waveforms back-propagated through a velocity model | No: hidden in the propagator | Partly: amplitude and radiation pattern consistency if elastic | Poorly: needs aperture and density |
| 4 Gradiometry, 6C, DAS | Spatial gradients, rotation, strain | Locally yes; it measures the slowness vector, i.e. the gradient of the travel-time field, without picks | Yes: direction and apparent velocity at one site | No: needs hardware not present |
| 5 Matched field, templates, back-projection | Cross-spectral or correlation match to modelled or empirical Green's functions | No: hidden in Green's function phase | Only if amplitudes are kept and trusted | Partly: empirical versions work, tied to known source regions |
| 6 Envelope, coda, energy | Envelopes, coda decay, energy | Envelope correlation: no (differential times). Energy and diffusion models: mostly yes | Yes: energy decay with distance; robust to overlap | Yes at high frequency, coarse resolution |
| 7 Machine learning | Raw waveforms | No: learned implicitly, with path fingerprints | Implicitly uses S-P, polarization, coda shape | Yes, but region-locked and outlier-prone |
| 8 Inverse-problem theory | Boundary (Cauchy) data, multi-frequency | Not separable | Framework, not a new observable | Theory assumes full boundary data |

## 2. Family by family

### 2.1 Polarization and particle motion

Observable and physics. The first P motion is rectilinear along the ray direction at the
receiver. The horizontal projection gives back azimuth, the vertical-to-horizontal ratio gives an
apparent incidence angle. Background: the apparent incidence angle differs from the true one
through the free-surface interaction and depends on the near-surface Vp/Vs, so incidence is a
weaker datum than azimuth.

Key references.

- Frohlich and Pulliam 1999 [M]: review of single-station location and a call for more work on
  it. Older than the window but the standard entry point.
- Agius and Galea 2011 [M], with the underlying thesis Agius 2007 [A]: automated single-station
  location for local and regional events in the Sicily Channel. Azimuth from three-component
  polarization, distance from S minus P. Reported accuracy: 70 percent of test events had both
  picks correct and 70 percent had azimuth within 30 degrees of bulletin values.
- Bose et al. 2017 [M]: probabilistic framework for single-station location, built for InSight.
  Background: it combines a distance pdf from differential phase times with an azimuth pdf from
  polarization.
- Zenhausern et al. 2022 [A]: time-frequency polarization analysis on P and S to obtain back
  azimuth for marsquakes, masking unpolarized signal; 24 events located, 16 of which had no
  previous back azimuth.
- Eisermann, Ziv and Wust-Bloch 2015 [M] and Netanel, Eisermann and Ziv 2021 [A]: real-time back
  azimuth for early warning. The 2021 paper locates off-network events from the P slowness vector
  of a small array plus S minus P, and intersects several arrays when available.

Geometry needed. One three-component station gives a direction. Two or more stations with
different azimuths to the source give an epicentre by intersection with no timing at all.

Accuracy. The only number read in this scan is the 30 degree figure above. Background: 10 to 20
degrees is the usual quoted range for good sites and impulsive P; this was not verified here.

Failure modes. Background: low signal-to-noise on horizontals, emergent P, lateral heterogeneity
and anisotropy that rotate the polarization away from the great circle, sensor misorientation,
and the 180 degree ambiguity when first-motion polarity is unclear. Cristiano, Schweitzer and
Meier 2013 [A] show at NORSAR that P polarization deviations are systematic, frequency dependent
and tied to dipping interfaces and anisotropy; systematic means they can be calibrated per station
and azimuth, which matters for the thesis.

Verdict. Truly independent of travel times. This is the cleanest non-kinematic datum available on
an existing three-component network.

### 2.2 Amplitude-based location

Observable and physics. Band-limited amplitude at each station, corrected for site amplification,
fitted to a decay law: geometric spreading times exponential attenuation.

Key references.

- Battaglia and Aki 2003 [A]: amplitudes corrected with coda site amplification factors; decay
  modelled as body or surface waves in a homogeneous medium; locates rockfalls, long-period events
  and eruption tremor at Piton de la Fournaise.
- Kumagai et al. 2010 [A]: the amplitude source location (ASL) method assuming isotropic S
  radiation; a band of about 5 to 12 Hz and Q of about 60 were found appropriate. Stated to work
  with few stations.
- Morioka, Kumagai and Maeda 2017 [A]: why the isotropic assumption works. Finite-difference
  simulations in heterogeneous structure show ASL is valid when the waveforms are multiply
  scattered and well described by a diffusion model, which reproduces the assumed
  amplitude-distance relation.
- Kumagai et al. 2019 [A]: ASL with depth-dependent scattering and attenuation at Nevado del
  Ruiz, Colombia; a two-layer attenuation model gives more accurate locations than a homogeneous
  one.
- Taisne et al. 2011 [M] and Caudron et al. 2018 [A]: inter-station amplitude ratios (SARA) to
  track migrating sources. Ratios cancel the source term. Caudron et al. find agreement with
  earthquake locations during the Bardarbunga intrusion, with exceptions, and state the method
  needs reliable site amplification and velocity information.
- Hardebeck and Shearer 2003 [M] and Funabiki and Miyazawa 2025 [A]: S/P amplitude ratios as a
  radiation-pattern constraint on focal mechanisms of small events. The 2025 paper obtains
  mechanisms from S/P ratios on a single-component DAS cable.

Geometry needed. Works with a handful of stations, provided site terms are calibrated and the
stations span a range of distances.

Accuracy. No general figure was read. Background: errors of one to a few kilometres are typical
inside a volcano network with apertures of 10 to 20 km; outside the network the amplitude decay
is nearly degenerate between distance and source size.

Failure modes. Site amplification errors map directly into location. Radiation pattern is not
isotropic at low frequency or short distance, so the isotropic model fails exactly where
waveform coherence is good. Two overlapping sources produce an amplitude pattern that is the sum
of two decays and is hard to separate without more stations.

Verdict. Truly independent of travel times. Two different uses must be kept apart: amplitude
decay with distance (needs incoherent, scattered high-frequency energy and calibrated site terms)
and radiation pattern (needs coherent, low-frequency amplitudes and is already inside any
moment-tensor fit).

### 2.3 Wave-equation and time-reversal imaging

Observable and physics. Recorded waveforms are reversed in time and injected at the receivers
into a numerical model; by reciprocity and time-reversal invariance the field refocuses at the
source. Variants differ in the imaging condition and in whether the misfit is minimised
explicitly (adjoint, full-waveform source inversion).

Key references.

- Larmat et al. 2006 [M] and Larmat et al. 2008 [A]: global time-reversal imaging with a 3D
  spectral-element propagator; the 2008 paper locates glacial earthquakes with a single
  time-reversed simulation instead of testing each candidate location.
- Gajewski and Tessmer 2005 [M], Artman, Podladtchikov and Witten 2010 [A]: reverse modelling
  for passive data. Artman et al. use elastic propagators, separate P and S in the image domain,
  use zero-lag auto- and cross-correlation imaging conditions, image the radiation pattern, and
  state that the result depends heavily on acquisition geometry.
- Nakata and Beroza 2016 [A]: geometric-mean imaging condition, cross-correlating the wavefields
  back-propagated from each receiver; reduces the search from four dimensions to three and
  suppresses side lobes.
- Kawakatsu and Montagner 2008 [M] and Kim, Liu and Tromp 2011 [M]: Background: these connect
  time reversal with moment-tensor inversion; time reversal is the adjoint operation, not the
  inverse, of the forward problem.
- Kremers et al. 2011 [A]: limits. Time-reversal imaging is useful only for effectively point
  sources; depth resolution is degraded by strong surface waves, incomplete surface coverage and
  smooth models that remove helpful scattering.
- Werner and Saenger 2018 [A]: systematic finite-difference study of station spacing, aperture,
  velocity-model error and noise; needs no event identification or onset picking.
- Douma and Snieder 2015 [A] and Bazargani and Snieder 2016 [A]: improved focusing by
  deconvolution, and optimised back-propagation signals for incomplete acquisition; the latter
  needs an approximate location but not the mechanism or origin time, and recovers double-couple
  radiation patterns better than plain time reversal.
- Sharan, Wang and Herrmann 2019 [A]: sparsity-promoting wave-equation source inversion, L1 in
  space and L2 in time; resolves closely spaced sources better than time reversal, recovers
  locations and source-time functions jointly without picks or a known event count.
- Sjogreen and Petersson 2014 [A]: full-waveform inversion for 11 point-source parameters (six
  moment-tensor components, location, origin time, a frequency parameter) with discrete adjoints;
  the problem is ill conditioned and needs Hessian scaling.
- Huang and Alkhalifah 2024 [A]: physics-informed neural network solving the frequency-domain
  wave equation with hard boundary constraints to image single and multiple sources from very
  sparse recordings, with fewer artefacts than time reversal.

Geometry needed. Aperture and density. Werner and Saenger quantify the requirement; a 4 to 8
station network with one-sided coverage is far from it.

Failure modes. Velocity-model error defocuses the image. Sparse receivers give strong artefacts,
because each receiver back-propagates an isochrone-like shell and too few shells intersect.
Surface waves dominate and blur depth.

Verdict. Hides travel times. The phase of the back-propagated field is the travel-time
information; a correct image needs the same kinematic accuracy as migration. What it adds, in the
elastic case, is consistency of amplitude, polarization and radiation pattern with the full
Green's tensor. A moment-tensor subspace detector built on full Green's functions already uses
that information, so this family adds little that is new to the thesis. Its useful ideas are the
imaging conditions (geometric mean, P-S cross-correlation) and the sparsity formulation of Sharan
et al., which is the closest published relative of a dictionary detector for overlapping events.

### 2.4 Gradiometry, rotational and strain seismology, DAS

Observable and physics. For a plane wave, the spatial gradient of displacement is proportional to
its time derivative times the slowness vector. Rotation rate and translational acceleration at one
point are therefore related through the local phase velocity and the propagation direction.

Key references.

- Langston 2007a, 2007b [M]: wave gradiometry. Background: horizontal slowness, azimuth, geometric
  spreading and radiation-pattern change are obtained from a small array a fraction of a
  wavelength across.
- Igel et al. 2007 [M]: ring-laser observations of earthquake-induced rotation.
- Sollberger et al. 2018 [M] and Schmelzbach et al. 2018 [A]: six-component (6C) polarization
  analysis. Schmelzbach et al. state that a single 6C station allows array-type processing:
  wavefield separation, isolation of S, local slowness estimation.
- Yuan et al. 2020 [A] and Yuan et al. 2021 [A]: portable rotation sensors validated on local
  earthquakes; back azimuth from combined rotation and translation at one site, used to track
  rupture and moving sources with S waves.
- Donner, Bernauer and Igel 2016 [A]: adding rotation to moment-tensor inversion. Half as many 6C
  stations match twice as many 3C stations; depth resolution improved by 175 percent and
  information gain by 53 percent in the synthetic test, mainly on Mxz, Myz and centroid depth.
- Lindsey et al. 2017 [A]: DAS arrays record earthquake wavefields; an L-shaped array gives back
  azimuth and slowness.
- van den Ende and Ampuero 2021 [A]: DAS beamforming performs poorly because strain is dominated
  by shallow scattered waves; converting strain rate to particle velocity helps.
- Li and Zhan 2018 [A]: DAS template matching finds events below the noise level; about 1500
  sensors give 90 percent of the detection performance.
- Ben-Zeev and Lior 2025 [A]: epicentre from beamforming on many short DAS subarrays plus phase
  time differences, resolving the ambiguity of a linear fibre.

Geometry needed. One 6C station, or one small dense array, or a fibre with changes of direction.

Verdict. Avoids inter-station travel times, but the slowness vector is by definition the gradient
of the travel-time field, so this is local kinematic information obtained without picks. It is
independent of network arrival times and of polarization (it also works for S and surface waves).
The Donner et al. result, that rotation resolves the depth-sensitive moment-tensor components, is
the most relevant to a moment-tensor detector. None of this is available on an existing sparse
three-component network. It belongs in the thesis as an outlook or a synthetic design study, not
as a method for the present data.

### 2.5 Matched field processing, empirical Green's functions, templates, back-projection

Observable and physics. The cross-spectral density matrix of the array data is compared with
replica vectors computed for trial sources (Bartlett: linear projection; MVDR: adaptive,
higher resolution, less robust). Empirical versions replace the modelled replica with data from
known events.

Key references.

- Schippkus and Hadziioannou 2022 [A]: review and method; numerically computed Green's functions
  for the full elastic wavefield replace analytical ones; applications to Southern California
  earthquakes and secondary microseisms; tutorial code.
- Harris and Kvaerna 2010 [M]: empirical matched field processing with seismic arrays.
  Background: replicas are calibrated on previous events, which absorbs unmodelled structure and
  gives resolution below the classical array limit, at the price of only covering regions with
  prior events.
- Gibbons and Ringdal 2006 [A]: array-stacked waveform correlation against templates; about 0.4
  magnitude units better detection than single-channel correlation.
- Yoon et al. 2015 [M]: FAST, a fingerprint and similarity-search detector. Background: it finds
  repeating waveforms without templates and does not locate.
- Ishii et al. 2005 [A] and Kiser and Ishii 2017 [M]: back-projection of P waves over a dense
  array.
- Vera Rodriguez, Sacchi and Gu 2012 [M]: sparse representation over a dictionary of Green's
  functions to recover origin time, hypocentre and moment tensor together. Background: this is
  the nearest published ancestor of a moment-tensor dictionary detector, and it inherits the grid.

Verdict. Hides travel times. A phase-only Bartlett processor is frequency-domain delay-and-sum;
back-projection is stacking along travel-time curves by construction. Amplitude-preserving matched
field with full Green's functions does add amplitude and radiation-pattern information, but only
if the Green's function amplitudes are trusted, which on real sparse networks they are not without
empirical calibration. Templates and empirical replicas avoid the velocity model entirely but can
only find events near known ones.

### 2.6 Envelope, coda and energy methods

Key references.

- Obara 2002 [M] and Wech and Creager 2008 [M]: tremor location by envelope cross-correlation.
  The enveloc documentation [A], which implements Wech and Creager, states the mechanism: a grid
  search that matches inter-station envelope correlation lags to predicted differential travel
  times. This is travel-time location with lags in place of picks.
- Mizuno and Ide 2019 [A]: envelope correlation recast as maximum likelihood with variance
  weighting; detects several simultaneous tremor sources as local maxima of the objective;
  horizontal 2 to 8 Hz envelopes from 313 stations; simultaneous sources are separated when about
  100 km apart; bootstrap errors above 2 km are discarded.
- Sawazaki and Enescu 2014 [A]: inversion of continuous high-frequency envelopes with radiative
  transfer theory for the space-time distribution of energy release during a mainshock and its
  early aftershocks, where events overlap and picks are impossible. The energy release follows
  an Omori-like decay and matches catalogue magnitudes; locations agree approximately with the
  aftershock distribution.
- Morioka, Kumagai and Maeda 2017 [A]: the diffusion model as the justification of
  amplitude-distance laws (see 2.2).

Verdict. Mixed. Envelope cross-correlation hides travel times. Energy-based envelope inversion
does not need picks, tolerates overlapping events by design, and uses a physically different
part of the signal (scattered energy and its decay), but its spatial resolution is coarse
("approximately" in Sawazaki and Enescu; 100 km separation in Mizuno and Ide) and it is not a
route to resolving a 5 to 10 km ambiguity. S minus P and other differential times avoid absolute
origin time, not travel times; they are useful because they are insensitive to origin-time
trade-off and to common-mode timing errors, not because they carry different physics.

### 2.7 Machine learning from waveforms

Key references.

- Perol, Gharbi and Denolle 2018 [A]: CNN on single-station waveforms; detection and
  classification into source regions; many more detections than the catalogue in Oklahoma.
- Mousavi and Beroza 2020 [A]: two Bayesian networks estimate epicentral distance, P travel time
  and back azimuth from one station; mean errors 7.3 km in epicentre, 0.4 s in origin time and
  6.7 km in depth.
- Kriegerowski et al. 2019 [M]: CNN on multistation full waveforms for clustered earthquakes.
- Zhang et al. 2020 [A]: fully convolutional network on 30 stations in Oklahoma; mean epicentre
  error 3.7 to 6.4 km for ML 1.5 and above, more than 11 km for smaller events; accuracy depends
  on the catalogue used as ground truth and on training coverage.
- van den Ende and Ampuero 2020 [A]: graph neural network that takes station coordinates as
  input; outperforms location-agnostic networks; works with a variable number of stations.
- Munchmeyer et al. 2021 [A]: transformer on a variable station set; epicentral error falls from
  7-22 km to 2-14 km within 25 s. The authors report heavy-tailed errors with outliers of hundreds
  of kilometres outside the training area, and state that location relies on region-specific
  scattering fingerprints, not on physics-based travel-time analysis, which limits transfer.
- Yang et al. 2021 [A] and Huang and Alkhalifah 2024 [A]: neural operators as fast surrogate
  wave solvers (2D acoustic in Yang et al.), and physics-informed networks for source imaging.

What physics they implicitly use. Background: single-station networks can only obtain distance
from S minus P and coda duration, and azimuth from polarization, so they are a learned version of
family 1 plus a path fingerprint. The fingerprint claim is not background: it is the authors' own
diagnosis in Munchmeyer et al. 2021.

Verdict. They do not avoid travel times; they learn them together with everything else and do
not say which datum drove the answer. Reported accuracy on sparse or single-station input (about
4 to 7 km at best, with heavy tails) is the same size as the mislocation the thesis is trying to
remove. Not a route to a mathematically grounded contribution, except as surrogate forward
solvers.

### 2.8 Inverse-problem formulations

Key references.

- El Badia and Ha-Duong 2000 [M] and El Badia and Nara 2011 [A]: algebraic identification of
  point sources from boundary data. The 2011 abstract states uniqueness and local stability for
  several point sources of the Helmholtz equation in an interior domain from Cauchy data at one
  wavenumber, and an algebraic algorithm for their number, positions and intensities.
  Background: the mechanism is the reciprocity gap functional tested against a family of
  solutions of the homogeneous equation.
- Bao et al. 2018 [A]: time-dependent elastodynamic sources that are a product of a spatial and a
  compactly supported temporal function; uniqueness of the spatial part from the field on a large
  sphere over a finite time; stability for the temporal part from one receiver; uniqueness with
  partial boundary data; Landweber and non-iterative schemes.
- Bao, Li and Zhao 2020 [A]: stability for elastic and electromagnetic inverse source problems;
  multi-frequency Dirichlet data overcome ill-posedness, with a stability estimate whose
  high-frequency tail shrinks as the band widens.
- Wang, Guo and Bousba 2021 preprint [A]: direct sampling for moment-tensor point sources of the
  Lame system; locations from maxima of an indicator function, then moment tensors from the
  intensities; multi-frequency displacement and traction on a surrounding boundary; sources must
  be well separated, and the separation fixes the lowest usable frequency; no forward solves.
- Kirsch 2002 [M]: MUSIC and the factorization method. Background: MUSIC locates point scatterers
  or sources from the noise subspace of a multistatic or multi-snapshot data matrix.
- Bredies and Pikkarainen 2013 [M], Candes and Fernandez-Granda 2014 [M], Duval and Peyre 2015
  [A], Pieper et al. 2020 [A]: sparse recovery in spaces of measures. Duval and Peyre prove that
  under a non-degenerate dual certificate and small noise the number, positions and amplitudes of
  spikes are recovered, and relate solutions on a grid to the continuous ones: the grid solution
  concentrates near the true positions rather than on them. Pieper et al. treat point sources of
  the Helmholtz equation in a bounded domain from few-frequency pointwise measurements with
  weights that grow near the observation points, solved by accelerated conditional gradient.
- Chi et al. 2011 [M]: sensitivity of compressed sensing to basis mismatch. Background: when the
  true parameter lies off the dictionary grid, the representation is no longer sparse and the
  energy leaks to neighbouring atoms.

Verdict. These do not separate kinematic from dynamic information; everything is in the wavefield
on the boundary. The theorems assume closed or large-aperture boundaries and many frequencies,
which a 4 to 8 station network does not provide, so they cannot be invoked as guarantees. Two
things transfer. First, the framing: with sparse sensors the point-source problem is not unique,
and what can be claimed is a resolution statement. Second, the off-grid results, which describe
the failure mode already observed in the thesis.

## 3. What is actually independent of arrival times

Strictly independent data that an existing three-component network records:

1. P polarization azimuth at each station (and, with more care, incidence).
2. Amplitude ratios between phases at one station (S/P), which carry radiation pattern.
3. Amplitude ratios between stations, which carry distance decay and radiation pattern, and need
   site terms.
4. High-frequency energy and coda decay.

Independent but needing new hardware: wavefield gradients, rotation, strain.

Not independent (kinematics hidden in a propagator or in learned weights): time reversal, RTM,
adjoint and full-waveform source inversion, matched field, back-projection, envelope
cross-correlation, templates, end-to-end neural location.

A useful geometric fact for the thesis. Background, simple trigonometry: travel-time and S minus P
information constrain mainly the source-station distance, while polarization constrains the
direction transverse to the ray. A 5 to 10 km error transverse to the path of a station 30 km
away changes the back azimuth by about 9.5 to 18 degrees; at 60 km by about 5 to 9.5 degrees. With
single-station azimuth errors of 10 to 30 degrees, one station does not resolve it, but two or
three reliable stations at different azimuths, or a station close to the source, can. Along the
ray the same 5 km changes amplitude by roughly 17 percent at 30 km under 1/r spreading, which is
well below a site-factor scatter of 0.5 in natural log; amplitude decay cannot resolve this
ambiguity without calibrated site terms.

## 4. Ranked shortlist for the thesis

Evidence taken as given: a moment-tensor subspace detector with timing tolerance; misses are
mostly events captured by atoms 5 to 10 km from the truth; networks of 4 to 8 three-component
stations; polarization reliable at only some stations on the real network. From REPORT_C8 in this
folder: per-station and per-phase site factors with log-normal spread 0.5, and catalogue amplitude
terms did not help.

### Rank 1. Timing-invariant polarization term with per-station reliability

Timing tolerance removes exactly the information that separates neighbouring atoms. P
polarization direction does not change when a trace is shifted in time, so it survives the
tolerance, and it constrains the direction that travel times constrain least. The proposal is a
directional likelihood per station (for instance a von Mises or projected-normal term on the
azimuth of the principal axis of the P covariance, as in Bose et al. 2017 and Zenhausern et al.
2022), weighted by a reliability estimated per station from catalogue events, with a station and
azimuth dependent bias term because the deviations are systematic (Cristiano et al. 2013). It is
used to re-rank atoms within the tolerance neighbourhood, not to detect. This fits the fact that
only some stations are reliable: unreliable stations get a flat likelihood and cost nothing. The
risk is stated in section 3: at 40 to 60 km the angular separation of atoms 5 km apart is below
single-station precision, so the gain must be shown per geometry, and it may be small. The test
is cheap and falsifiable: on the existing misses, does the true atom have higher polarization
likelihood than the capturing atom more often than chance?

### Rank 2. Off-grid, measure-valued refinement of the dictionary

The observed failure, the right event on the wrong nearby atom, is what basis mismatch and the
grid results of Duval and Peyre predict for any gridded sparse method. Before adding new physics
it should be settled how much of the 5 to 10 km error is a grid effect and how much is a true
ambiguity of the sparse network. The direction is a continuous-location refinement step after
atom selection (local optimisation of location jointly with the moment tensor, as in conditional
gradient methods with sliding, Pieper et al. 2020; or the adjoint gradient of Sjogreen and
Petersson 2014), together with a resolution statement in place of a uniqueness claim. This has
the strongest mathematical content of the five and connects the thesis to an existing theory
(Candes and Fernandez-Granda; Duval and Peyre; Bredies and Pikkarainen). It is not a
travel-time-free method. If the ambiguity is geometric, refinement will not remove it and Rank 1
or Rank 3 is required; that outcome is itself a result.

### Rank 3. Phase amplitude ratios (S/P) as a robust radiation-pattern constraint

S/P ratios at one station cancel the source size and part of the path, and discriminate between
atoms whose take-off angles differ (Hardebeck and Shearer 2003; Funabiki and Miyazawa 2025).
A full-waveform moment-tensor fit already contains this information in principle, but a
timing-tolerant subspace match can dilute it, and REPORT_C8 found that a robust (Huber) treatment
of station and phase residual norms was the only principled variant that helped. The direction
is to make the ratio constraint explicit and robust. The limit is known from C8: site factors
were drawn per station and per phase, so S and P site terms do not cancel in the ratio; the gain
depends on whether real site terms are correlated between phases, which should be measured on
real data first. Inter-station amplitude decay (ASL proper) is ranked lower still for this
problem, for the reason computed in section 3.

### Rank 4. High-frequency envelope and energy channel for overlapping events

When events overlap, coherent waveform matching degrades and energy methods do not (Sawazaki and
Enescu 2014; Mizuno and Ide 2019; diffusion-model justification in Morioka et al. 2017). A second,
incoherent channel (band-limited envelopes fitted with a radiative-transfer or diffusion kernel)
could give the number of events and coarse energy centroids in a window, to be used as a prior or
a model-order check for the coherent detector. This addresses detection and counting under
overlap, not the 5 to 10 km mislocation; published resolution is far coarser than that. It also
needs site amplification from coda, which the Kumagai and Battaglia papers treat as a
prerequisite.

### Rank 5. Elastic time-reversal or geometric-mean focusing as a post-detection discriminator

Back-propagating the window through the same Green's functions and applying a multiplicative
imaging condition (Nakata and Beroza 2016) or an optimised back-propagation (Bazargani and
Snieder 2016) gives a sharper function of position than a linear stack, and the sparse
wave-equation formulation of Sharan et al. 2019 is the nearest published competitor that should
be cited and, if feasible, compared. It ranks last because it hides travel times, adds no datum
beyond what the Green's function atoms already encode, and its documented requirements on
aperture and density (Werner and Saenger 2018; Kremers et al. 2011) are not met by 4 to 8
stations.

Not shortlisted. Gradiometry, 6C and DAS: genuinely independent information and a strong
published case for moment tensors (Donner et al. 2016), but the hardware is absent; suitable as an
outlook section or a synthetic network-design study. End-to-end machine learning: published
single-station and sparse-network accuracy is the same size as the error to be removed, with
heavy tails and no transfer across regions (Munchmeyer et al. 2021).

## 5. Things I could not establish

- A verified, general accuracy figure for P polarization back azimuth on local earthquakes. The
  only number read is 70 percent within 30 degrees (Agius 2007).
- A verified accuracy figure for ASL in kilometres.
- Whether any paper combines a moment-tensor dictionary detector with an explicit polarization
  likelihood. I found none, but the search was limited by rate limits and by search results that
  returned titles only; absence here is not evidence of novelty.
- Coverage of 2025-2026 is thin: Funabiki and Miyazawa 2025 and Ben-Zeev and Lior 2025 are the
  only items from those years.

## 6. References

Tags: [A] abstract or page read; [M] metadata only. "Opened at" gives the page actually read.

Polarization and single station

1. Frohlich, C., and Pulliam, J. (1999). Single-station location of seismic events: a review and
   a plea for more research. Physics of the Earth and Planetary Interiors 113, 277-291.
   doi:10.1016/S0031-9201(99)00055-2. [M] Opened at api.crossref.org.
2. Agius, M. R., and Galea, P. (2011). A single-station automated earthquake location system at
   Wied Dalam station, Malta. Seismological Research Letters 82, 545-559.
   doi:10.1785/gssrl.82.4.545. [M] Opened at api.crossref.org.
3. Agius, M. R. (2007). Automatic earthquake detection and localisation from a three-component
   single-station. Master's dissertation, University of Malta.
   https://www.um.edu.mt/library/oar/handle/123456789/79083 [A]
4. Bose, M., Clinton, J. F., Ceylan, S., Euchner, F., van Driel, M., Khan, A., Giardini, D.,
   Lognonne, P., and Banerdt, W. B. (2017). A probabilistic framework for single-station location
   of seismicity on Earth and Mars. Physics of the Earth and Planetary Interiors 262, 48-65.
   doi:10.1016/j.pepi.2016.11.003. [M] Opened at api.crossref.org.
5. Zenhausern, G., Stahler, S. C., Clinton, J. F., Giardini, D., Ceylan, S., and Garcia, R. F.
   (2022). Low-frequency marsquakes and where to find them: back azimuth determination using a
   polarization analysis approach. Bulletin of the Seismological Society of America 112,
   1787-1805. doi:10.1785/0120220019. [A] Opened at api.crossref.org.
6. Eisermann, A. S., Ziv, A., and Wust-Bloch, G. H. (2015). Real-time back azimuth for earthquake
   early warning. Bulletin of the Seismological Society of America 105, 2274-2285.
   doi:10.1785/0120140298. [M] Opened at api.crossref.org.
7. Netanel, M., Eisermann, A. S., and Ziv, A. (2021). Off-network earthquake location by
   earthquake early warning systems: methodology and validation. Bulletin of the Seismological
   Society of America 111, 3090-3102. doi:10.1785/0120210006. [A] Opened at cris.iucc.ac.il.
8. Cristiano, L., Schweitzer, J., and Meier, T. (2013). P-wave polarization and FK analysis of
   NORSAR array data to investigate local anisotropy and lateral heterogeneity. NORSAR Scientific
   Report 1-2013, 57-63. No DOI.
   https://www.norsar.no/getfile.php/1318446-1681728350/norsar.no/R-D/Publications/4979_Cristiano_etal2013.pdf [A]

Amplitude

9. Battaglia, J., and Aki, K. (2003). Location of seismic events and eruptive fissures on the
   Piton de la Fournaise volcano using seismic amplitudes. Journal of Geophysical Research 108
   (B8), 2364. doi:10.1029/2002JB002193. [A] Opened at pubs.usgs.gov/publication/70025168 and
   api.crossref.org. The article number 2364 was not shown on the pages read.
10. Kumagai, H., Nakano, M., Maeda, T., Yepes, H., Palacios, P., Ruiz, M., Arrais, S., Vaca, M.,
    Molina, I., and Yamashima, T. (2010). Broadband seismic monitoring of active volcanoes using
    deterministic and stochastic approaches. Journal of Geophysical Research 115, B8.
    doi:10.1029/2009JB006889. [A] Opened at api.crossref.org.
11. Kumagai, H., Lacson, R., Maeda, Y., Figueroa, M. S., Yamashina, T., Ruiz, M., Palacios, P.,
    Ortiz, H., and Yepes, H. (2013). Source amplitudes of volcano-seismic signals determined by
    the amplitude source location method as a quantitative measure of event size. Journal of
    Volcanology and Geothermal Research 257, 57-71. doi:10.1016/j.jvolgeores.2013.03.002. [M]
12. Morioka, H., Kumagai, H., and Maeda, T. (2017). Theoretical basis of the amplitude source
    location method for volcano-seismic signals. Journal of Geophysical Research: Solid Earth
    122, 6538-6551. doi:10.1002/2017JB013997. [A] Opened at api.crossref.org.
13. Kumagai, H., Londono, J. M., Maeda, Y., and Acevedo Rivas, A. E. (2019). Amplitude source
    location method with depth-dependent scattering and attenuation structures: application at
    Nevado del Ruiz volcano, Colombia. Journal of Geophysical Research: Solid Earth 124,
    11585-11600. doi:10.1029/2019JB018156. [A] Opened at api.crossref.org.
14. Taisne, B., Brenguier, F., Shapiro, N. M., and Ferrazzini, V. (2011). Imaging the dynamics of
    magma propagation using radiated seismic intensity. Geophysical Research Letters 38.
    doi:10.1029/2010GL046068. [M] Opened at api.crossref.org.
15. Caudron, C., White, R. S., Green, R. G., Woods, J., Agustsdottir, T., Donaldson, C.,
    Greenfield, T., Rivalta, E., and Brandsdottir, B. (2018). Seismic amplitude ratio analysis of
    the 2014-2015 Bardarbunga-Holuhraun dike propagation and eruption. Journal of Geophysical
    Research: Solid Earth 123, 264-276. doi:10.1002/2017JB014660. [A] Opened at opinvisindi.is.
16. Hardebeck, J. L., and Shearer, P. M. (2003). Using S/P amplitude ratios to constrain the focal
    mechanisms of small earthquakes. Bulletin of the Seismological Society of America 93,
    2434-2444. doi:10.1785/0120020236. [M] Opened at api.crossref.org (the record listed only the
    first author; the second author is from memory).
17. Funabiki, Y., and Miyazawa, M. (2025). Estimating focal mechanism of small earthquakes using
    S/P amplitude ratios of distributed acoustic sensing records. Geophysical Research Letters
    52, e2024GL113963. doi:10.1029/2024GL113963. [A] Opened at api.crossref.org.

Time reversal, wave-equation imaging, waveform source inversion

18. Larmat, C., Montagner, J.-P., Fink, M., Capdeville, Y., Tourin, A., and Clevede, E. (2006).
    Time-reversal imaging of seismic sources and application to the great Sumatra earthquake.
    Geophysical Research Letters 33. doi:10.1029/2006GL026336. [M] Opened at
    api.semanticscholar.org.
19. Larmat, C., Tromp, J., Liu, Q., and Montagner, J.-P. (2008). Time reversal location of glacial
    earthquakes. Journal of Geophysical Research 113, B9. doi:10.1029/2008JB005607. [A]
20. Gajewski, D., and Tessmer, E. (2005). Reverse modelling for seismic event characterization.
    Geophysical Journal International 163, 276-284. doi:10.1111/j.1365-246X.2005.02732.x. [M]
21. Artman, B., Podladtchikov, I., and Witten, B. (2010). Source location using time-reverse
    imaging. Geophysical Prospecting 58, 861-873. doi:10.1111/j.1365-2478.2010.00911.x. [A]
22. Nakata, N., and Beroza, G. C. (2016). Reverse time migration for microseismic sources using
    the geometric mean as an imaging condition. Geophysics 81, KS51-KS60.
    doi:10.1190/geo2015-0278.1. [A]
23. Kawakatsu, H., and Montagner, J.-P. (2008). Time-reversal seismic-source imaging and
    moment-tensor inversion. Geophysical Journal International 175, 686-688.
    doi:10.1111/j.1365-246X.2008.03926.x. [M]
24. Kim, Y., Liu, Q., and Tromp, J. (2011). Adjoint centroid-moment tensor inversions.
    Geophysical Journal International. doi:10.1111/j.1365-246X.2011.05027.x. [M] Opened at
    api.semanticscholar.org.
25. Kremers, S., Fichtner, A., Brietzke, G. B., Igel, H., Larmat, C., Huang, L., and Kaser, M.
    (2011). Exploring the potentials and limitations of the time-reversal imaging of finite
    seismic sources. Solid Earth 2, 95-105. doi:10.5194/se-2-95-2011. [A] Opened at
    se.copernicus.org.
26. Werner, C., and Saenger, E. H. (2018). Obtaining reliable source locations with time reverse
    imaging: limits to array design, velocity models and signal-to-noise ratios. Solid Earth 9,
    1487-1505. doi:10.5194/se-9-1487-2018. [A]
27. Douma, J., and Snieder, R. (2015). Focusing of elastic waves for microseismic imaging.
    Geophysical Journal International 200, 390-401. doi:10.1093/gji/ggu398. [A] Opened at
    academic.oup.com.
28. Bazargani, F., and Snieder, R. (2016). Optimal source imaging in elastic media. Geophysical
    Journal International 204, 1134-1147. doi:10.1093/gji/ggv494. [A] Opened at academic.oup.com.
29. Sharan, S., Wang, R., and Herrmann, F. J. (2019). Fast sparsity-promoting microseismic source
    estimation. Geophysical Journal International 216, 164-181. doi:10.1093/gji/ggy415. [A]
30. Sjogreen, B., and Petersson, N. A. (2014). Source estimation by full wave form inversion.
    Journal of Scientific Computing 59, 247-276. doi:10.1007/s10915-013-9760-6. [A] Opened at
    link.springer.com.
31. Huang, X., and Alkhalifah, T. (2024). Microseismic source imaging using physics-informed
    neural networks with hard constraints. IEEE Transactions on Geoscience and Remote Sensing.
    doi:10.1109/TGRS.2024.3366449. arXiv:2304.04315. [A] Opened at arxiv.org/abs/2304.04315.

Gradiometry, rotation, DAS

32. Langston, C. A. (2007a). Wave gradiometry in two dimensions. Bulletin of the Seismological
    Society of America 97, 401-416. doi:10.1785/0120060138. [M]
33. Langston, C. A. (2007b). Wave gradiometry in the time domain. Bulletin of the Seismological
    Society of America 97, 926-933. doi:10.1785/0120060152. [M]
34. Igel, H., Cochard, A., Wassermann, J., Flaws, A., Schreiber, U., Velikoseltsev, A., and Pham
    Dinh, N. (2007). Broad-band observations of earthquake-induced rotational ground motions.
    Geophysical Journal International. doi:10.1111/j.1365-246X.2006.03146.x. [M]
35. Sollberger, D., Greenhalgh, S. A., Schmelzbach, C., Van Renterghem, C., and Robertsson,
    J. O. A. (2018). 6-C polarization analysis using point measurements of translational and
    rotational ground-motion: theory and applications. Geophysical Journal International 213,
    77-97. doi:10.1093/gji/ggx542. [M] Opened at api.openalex.org.
36. Schmelzbach, C., Donner, S., Igel, H., Sollberger, D., Taufiqurrahman, T., Bernauer, F.,
    Hausler, M., Van Renterghem, C., Wassermann, J., and Robertsson, J. (2018). Advances in 6C
    seismology: applications of combined translational and rotational motion measurements in
    global and exploration seismology. Geophysics 83, WC53-WC69. doi:10.1190/geo2017-0492.1. [A]
37. Yuan, S., Simonelli, A., Lin, C.-J., Bernauer, F., Donner, S., Braun, T., Wassermann, J., and
    Igel, H. (2020). Six degree-of-freedom broadband ground-motion observations with portable
    sensors: validation, local earthquakes, and signal processing. Bulletin of the Seismological
    Society of America 110, 953-969. doi:10.1785/0120190277. [A]
38. Yuan, S., Gessele, K., Gabriel, A.-A., May, D. A., Wassermann, J., and Igel, H. (2021).
    Seismic source tracking with six degree-of-freedom ground motion observations. Journal of
    Geophysical Research: Solid Earth 126, e2020JB021112. doi:10.1029/2020JB021112. [A]
39. Donner, S., Bernauer, M., and Igel, H. (2016). Inversion for seismic moment tensors combining
    translational and rotational ground motions. Geophysical Journal International 207, 562-570.
    doi:10.1093/gji/ggw298. [A] Opened at academic.oup.com.
40. Lindsey, N. J., Martin, E. R., Dreger, D. S., and others (2017). Fiber-optic network
    observations of earthquake wavefields. Geophysical Research Letters 44 (23).
    doi:10.1002/2017GL075722. [A]
41. van den Ende, M. P. A., and Ampuero, J.-P. (2021). Evaluating seismic beamforming
    capabilities of distributed acoustic sensing arrays. Solid Earth 12, 915-934.
    doi:10.5194/se-12-915-2021. [A]
42. Li, Z., and Zhan, Z. (2018). Pushing the limit of earthquake detection with distributed
    acoustic sensing and template matching: a case study at the Brady geothermal field.
    Geophysical Journal International 215, 1583-1593. doi:10.1093/gji/ggy359. [A]
43. Ben-Zeev, S., and Lior, I. (2025). Earthquake location with distributed acoustic sensing
    subarray beamforming with implications for earthquake early warning. Seismological Research
    Letters 96, 2453-2462. doi:10.1785/0220240452. [A] Opened at cris.huji.ac.il.

Matched field, templates, back-projection, dictionaries

44. Schippkus, S., and Hadziioannou, C. (2022). Matched field processing accounting for complex
    Earth structure: method and review. Geophysical Journal International 231, 1268-1282.
    doi:10.1093/gji/ggac240. [A] Opened at api.openalex.org.
45. Harris, D. B., and Kvaerna, T. (2010). Superresolution with seismic arrays using empirical
    matched field processing. Geophysical Journal International 182, 1455-1477.
    doi:10.1111/j.1365-246X.2010.04684.x. [M]
46. Gibbons, S. J., and Ringdal, F. (2006). The detection of low magnitude seismic events using
    array-based waveform correlation. Geophysical Journal International 165, 149-166.
    doi:10.1111/j.1365-246X.2006.02865.x. [A]
47. Yoon, C. E., O'Reilly, O., Bergen, K. J., and Beroza, G. C. (2015). Earthquake detection
    through computationally efficient similarity search. Science Advances 1, e1501057.
    doi:10.1126/sciadv.1501057. [M]
48. Ishii, M., Shearer, P. M., Houston, H., and Vidale, J. E. (2005). Extent, duration and speed
    of the 2004 Sumatra-Andaman earthquake imaged by the Hi-Net array. Nature 435, 933-936.
    doi:10.1038/nature03675. [A]
49. Kiser, E., and Ishii, M. (2017). Back-projection imaging of earthquakes. Annual Review of
    Earth and Planetary Sciences 45, 271-299. doi:10.1146/annurev-earth-063016-015801. [M]
50. Vera Rodriguez, I., Sacchi, M., and Gu, Y. J. (2012). Simultaneous recovery of origin time,
    hypocentre location and seismic moment tensor using sparse representation theory. Geophysical
    Journal International 188, 1188-1202. doi:10.1111/j.1365-246X.2011.05323.x. [M]

Envelope, coda, energy

51. Obara, K. (2002). Nonvolcanic deep tremor associated with subduction in southwest Japan.
    Science. doi:10.1126/science.1070378. [M]
52. Wech, A. G., and Creager, K. C. (2008). Automated detection and location of Cascadia tremor.
    Geophysical Research Letters. doi:10.1029/2008GL035458. [M] Method description read in the
    enveloc documentation, https://enveloc.readthedocs.io/en/stable/intro.html [A]
53. Mizuno, N., and Ide, S. (2019). Development of a modified envelope correlation method based
    on maximum-likelihood method and application to detecting and locating deep tectonic tremors
    in western Japan. Earth, Planets and Space 71, 40. doi:10.1186/s40623-019-1022-x. [A]
54. Sawazaki, K., and Enescu, B. (2014). Imaging the high-frequency energy radiation process of a
    main shock and its early aftershock sequence: the case of the 2008 Iwate-Miyagi Nairiku
    earthquake, Japan. Journal of Geophysical Research: Solid Earth 119, 4729-4746.
    doi:10.1002/2013JB010539. [A]

Machine learning

55. Perol, T., Gharbi, M., and Denolle, M. (2018). Convolutional neural network for earthquake
    detection and location. Science Advances. doi:10.1126/sciadv.1700578. [A] Opened at
    api.semanticscholar.org (which lists the year as 2017).
56. Mousavi, S. M., and Beroza, G. C. (2020). Bayesian-deep-learning estimation of earthquake
    location from single-station observations. IEEE Transactions on Geoscience and Remote
    Sensing. doi:10.1109/TGRS.2020.2988770. arXiv:1912.01144. [A]
57. Kriegerowski, M., Petersen, G. M., Vasyura-Bathke, H., and Ohrnberger, M. (2019). A deep
    convolutional neural network for localization of clustered earthquakes based on multistation
    full waveforms. Seismological Research Letters 90, 510-516. doi:10.1785/0220180320. [M]
58. Zhang, X., and others (2020). Locating induced earthquakes with a network of seismic stations
    in Oklahoma via a deep learning method. Scientific Reports 10, 1941.
    doi:10.1038/s41598-020-58908-5. [A] Opened at nature.com; the page summary gave "Zhang et
    al." without the full author list.
59. van den Ende, M. P. A., and Ampuero, J.-P. (2020). Automated seismic source characterization
    using deep graph neural networks. Geophysical Research Letters. doi:10.1029/2020GL088690. [A]
60. Munchmeyer, J., Bindi, D., Leser, U., and Tilmann, F. (2021). Earthquake magnitude and
    location estimation from real time seismic waveforms with a transformer network. Geophysical
    Journal International 226, 1086-1104. doi:10.1093/gji/ggab139. [A] Opened at academic.oup.com
    (the page summary gave "Munchmeyer et al."; co-authors are from memory).
61. Yang, Y., Gao, A. F., Castellanos, J. C., Ross, Z. E., Azizzadenesheli, K., and Clayton, R. W.
    (2021). Seismic wave propagation and inversion with neural operators. arXiv:2108.05421.
    https://arxiv.org/abs/2108.05421 [A] The journal version was not opened.

Inverse-problem theory and sparse recovery

62. El Badia, A., and Ha-Duong, T. (2000). An inverse source problem in potential analysis.
    Inverse Problems 16, 651-663. doi:10.1088/0266-5611/16/3/308. [M]
63. El Badia, A., and Nara, T. (2011). An inverse source problem for Helmholtz's equation from the
    Cauchy data with a single wave number. Inverse Problems 27, 105001.
    doi:10.1088/0266-5611/27/10/105001. [A] Opened at doi.org and api.semanticscholar.org.
64. Bao, G., Hu, G., Kian, Y., and Yin, T. (2018). Inverse source problems in elastodynamics.
    Inverse Problems 34, 045009. doi:10.1088/1361-6420/aaaf7e. arXiv:1708.04684. [A] Opened at
    arxiv.org; volume and article number are from memory.
65. Bao, G., Li, P., and Zhao, Y. (2020). Stability for the inverse source problems in elastic and
    electromagnetic waves. Journal de Mathematiques Pures et Appliquees 134, 122-178.
    doi:10.1016/j.matpur.2019.06.006. [A] Opened at math.purdue.edu (author copy).
66. Wang, X., Guo, Y., and Bousba, S. (2021). Direct imaging for the moment tensor point sources
    of elastic waves. arXiv:2101.10557. https://arxiv.org/abs/2101.10557 [A] The journal version,
    if any, was not opened.
67. Kirsch, A. (2002). The MUSIC-algorithm and the factorization method in inverse scattering
    theory for inhomogeneous media. Inverse Problems 18, 1025-1040.
    doi:10.1088/0266-5611/18/4/306. [M]
68. Bredies, K., and Pikkarainen, H. K. (2013). Inverse problems in spaces of measures. ESAIM:
    Control, Optimisation and Calculus of Variations 19, 190-218. doi:10.1051/cocv/2011205. [M]
69. Candes, E. J., and Fernandez-Granda, C. (2014). Towards a mathematical theory of
    super-resolution. Communications on Pure and Applied Mathematics 67, 906-956.
    doi:10.1002/cpa.21455. [M]
70. Duval, V., and Peyre, G. (2015). Exact support recovery for sparse spikes deconvolution.
    Foundations of Computational Mathematics 15, 1315-1355. doi:10.1007/s10208-014-9228-6. [A]
71. Pieper, K., Tang, B. Q., Trautmann, P., and Walter, D. (2020). Inverse point source location
    with the Helmholtz equation on a bounded domain. Computational Optimization and Applications
    77, 213-249. doi:10.1007/s10589-020-00205-y. [A]
72. Chi, Y., Scharf, L. L., Pezeshki, A., and Calderbank, A. R. (2011). Sensitivity to basis
    mismatch in compressed sensing. IEEE Transactions on Signal Processing 59, 2182-2195.
    doi:10.1109/TSP.2011.2112650. [M]
