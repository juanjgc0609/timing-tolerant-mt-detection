# Appendix AB. Phase 5 Theory: Identifiability of Overlapping Events, Greedy Recovery with a Perturbed Dictionary, an Adaptive Robust Statistic and the Components of the New Model

Date: 18 October 2026. Written against `main` at `a936cd6`, after the experiments, and integrated into the repository
on 19 October 2026 as the companion of Section 4 of the thesis (`docs/thesis/Thesis.md`); this file continues the
evidence volume `Supplementary_Evidence_Appendices_G_to_AA.md`. The scripts and their stored outputs are in
`experiments/exploratory/THEORY-5/`; they were re-run in a clean Linux environment before integration (see that
folder's `card.md`). Parts A to C bear on the thesis; Part D is theory for the doctoral line
(`docs/decisions/ADR-0006-two-products.md`) and no claim of the thesis rests on it.

Taxonomy, the same as in Appendix AA: a **Proposition**, **Lemma**, **Theorem** or **Corollary** is proved here
completely under the hypotheses it lists; a **Conjecture** is not proved; an **Approximation** is a formula used
without an error bound. Every figure comes from a script in `experiments/exploratory/THEORY-5/` with its output stored
next to it (`experiments/exploratory/THEORY-5/*_out.txt`); the scripts are self-contained
(`experiments/exploratory/THEORY-5/lib5.py` replicates the repository's conventions: 2 Hz Ricker, 16 Hz,
$V_S=3$ km/s, $V_P=1.73V_S$, 8-station `spread` network, radiation operator of `code/task1x_mtglrt.py`).

---

## AB.0 Summary

| # | Result | Status | Key figure |
|:--|:--|:--|:--|
| AB.A.1 | Le Cam bound: no test distinguishes "one" from "two events" with total error smaller than $2\Phi(-\delta/2)$ | Theorem | error 0.1 requires $\delta\ge3.29$ |
| AB.A.2 | Two coincident events are never identifiable, whatever the mechanisms | Proposition | — |
| AB.A.3 | The margin against the one-event class with **free per-arrival** times and amplitudes separates by station and is invariant to common statics | Theorem | — |
| AB.A.4–AB.A.5 | With different mechanisms the margin grows **linearly** with the separation; with the same mechanism, **quadratically** | Theorem | measured slopes 1.00 and 2.00 |
| AB.A.6 | Network form: $\delta^2\approx\omega^2\,\Delta\boldsymbol\theta^{\mathsf T}F_{\rm pol}\Delta\boldsymbol\theta$, observability weighted by polarization diversity | Corollary | — |
| AB.A.7 | Minimum resolvable separation: $\propto\mathrm{SNR}^{-1}$ with mechanism diversity, $\propto\mathrm{SNR}^{-1/2}$ without it | Corollary | 0.034 s and 0.098 s at SNR 12.4 |
| AB.B.1 | Exact condition for the energy-peak reading (D1) to read the correct arrival | Lemma | interference $<0.144\,\lVert\mathbf g\rVert$ |
| AB.B.2–AB.B.4 | Computable certificate of correct selection for the deployed algorithm with statics and overlap | Lemma and Theorem | correct in 100 % of cases, but it never activates: actual margin ≈15 % |
| AB.B.5 | Recovery of $K$ events by induction | Theorem | — |
| AB.B.6 | The deployed deflation leaves up to a fixed fraction of each arrival because of the sub-sample offset | Proposition | median residual 16 % (up to 36 %) of the amplitude |
| AB.C.1 | **Correction of AA.9.5(b)**: the Lipschitz constant is 1, not $\kappa=2.146$ | Proposition | bound 28.50 → 14.59 (empirical 7.12) |
| AB.C.2 | **Conditional** tail valid for any amplitude law (heavy tails, non-stationary, isotropic transients) | Theorem | median threshold 9.75 at $p=10^{-4}$ |
| AB.C.3 | **Impossibility**: influence bounded in the absolute sense implies saturation | Proposition | — |
| AB.C.4 | ACAP (cap relative to the block median): scale equivariant (does not saturate), bounded relative influence, zero cost if max/median $\le c_1$ | Proposition | — |
| AB.D.1 | Contamination likelihood: influence of a transient $\propto\log(1/f_{\rm background})$; the evidence for an event grows as $\nu\log a$ (tail index of the background) | Proposition | measured slope 3.00 with $\nu=3$ |
| AB.D.2 | The mixture version over the tensor is an e-value: universal threshold $\log(1/p)$ without calibration | Theorem | — |
| AB.D.3 | *Time-slides*: exact p-value for the global null "no coherent source", with no floor; does not control false detections caused by real events | Proposition | — |
| AB.D.4 | Joint decomposition: $\sigma_{\min}(\Gamma)\propto\Delta\tau$, so certificates of type AA.6(a) require SNR $\propto\Delta\tau^{-3}$ against the necessary $\Delta\tau^{-1}$ | Proposition | measured slope 2.00 of $1-\cos\theta_1$ |

Four messages for the plan:

1. **The information limit is not the bottleneck.** For a pair of events at the level of the Task 2M threshold
   (SNR 12.4 per event), the minimum separation resolvable by *any* method is 0.03 to 0.1 s in origin time
   or 0.5 km in plan view. The detector works with a tolerance of 0.25 s and a hit rule of 5 km and 1 s. The
   target regime of the atlas is set by model error and by the algorithm, not by identifiability.
2. **Timing tolerance costs little identifiability near overlap.** The margin against one-event explanations
   with *arbitrary* per-arrival times is 84 % of the exact margin at small separations, and it does not
   depend on the station statics (Theorem AB.A.3). What does cost is the D1 reading rule and the deflation on the
   sample grid (Lemma AB.B.1, Proposition AB.B.6).
3. **Bounded influence and non-saturation are compatible only if the bound is relative to the background** (Proposition
   AB.C.3). This gives mathematical form to recommendation H5 of the model design note (state of the art and design of
   the new model; kept in the project's shared files, not yet in the repository) and to its contamination
   likelihood (§4.4 of that note), and ACAP is the simplest radial case that satisfies it.
4. **For the new model**, the contamination likelihood does not saturate and its mixture version is an e-value with a
   universal threshold, but only if the background is modelled with polynomial tails (Proposition AB.D.1, Theorem AB.D.2); the
   *time-slides* calibrate the global null exactly and not the false detections that a real event causes (Proposition AB.D.3); and
   no coherence certificate can guarantee the separation of overlapping events near the information
   limit, by a factor $\Delta\tau^{-2}$ in the SNR (Proposition AB.D.4).

---

## Common notation

The data are whitened and sampled with step $\Delta=1/16$ s: $\mathbf d=\sum_k G(\theta_k)\mathbf m_k+\mathbf n$,
$\mathbf n\sim\mathcal N(0,I)$, with the arrival structure (AA.4): at station $i$ and branch $b\in\{P,S\}$, the
source $k$ contributes $\mathbf g^k_{ib}\,\bar s(\cdot-t^k_{ib})$, with $\mathbf g^k_{ib}=A_{ib}(\mathbf x_k)\mathbf m_k\in\mathbb R^3$ and
$\bar s$ a unit-energy pulse. The autocorrelation $a(t)=\langle\bar s,\bar s(\cdot-t)\rangle$ is even, $a(0)=1$,
$|a(t)|<1$ for $t\ne0$ and $a(t)\to0$. We write

$$\omega^2=\lVert\dot{\bar s}\rVert^2=-a''(0),\qquad \psi=\lVert\ddot{\bar s}\rVert^2=a^{(4)}(0).$$

For the Ricker with peak $f_c$ (power spectrum $\propto f^4e^{-2f^2/f_c^2}$) the autocorrelation has the closed form

$$a(t)=\Big(1-2u+\tfrac{u^2}{3}\Big)e^{-u/2},\quad u=(\pi f_ct)^2,\qquad \omega^2=5\pi^2f_c^2,\qquad \psi=35\pi^4f_c^4,$$

which is obtained by differentiating $\mathbb E\cos(2\pi tF)=e^{-\pi^2f_c^2t^2/2}$ four times for $F\sim\mathcal N(0,f_c^2/4)$.
With $f_c=2$ Hz: $\omega^2=197.4\ \mathrm s^{-2}$, $\psi-\omega^4=10\pi^4f_c^4=15{,}585\ \mathrm s^{-4}$, $a(0.5\,\mathrm s)=0.099$
(the 0.10 of AA.1), the deepest negative lobe is $a^-:=-\min a=0.618$ at $t=0.216$ s, and $a(\Delta/2)=0.906$. The
16 Hz sampling does not alter these inner products: by the Poisson summation formula the error is of the order of
$|\hat s(16-6\ \mathrm{Hz})|^2\sim e^{-50}$ (verified against the numerical autocorrelation to $10^{-10}$ in
`lib5.py`).

---

# Part A. Identifiability of overlapping events

## AB.A.0 Setting

Let $\xi=(\theta_1,\mathbf m_1,\theta_2,\mathbf m_2)$ be a pair of events and $\boldsymbol\mu_\xi=G(\theta_1)\mathbf m_1+G(\theta_2)\mathbf m_2$ its
signal. Two classes of one-event explanations:

- the **physical** one, $\mathcal M_1=\{G(\theta)\mathbf m:\ \theta\in\Theta,\ \mathbf m\in\mathbb R^6\}$;
- the **timing-tolerant** one, $\mathcal M_1^{\rm TT}=\{\sum_{i}\sum_{b}\mathbf e_i\otimes\mathbf h_{ib}\otimes\bar s(\cdot-t_{ib}):\ \mathbf h_{ib}\in\mathbb R^3,\ t_{ib}\in\mathbb R\}$,
  one pulse per station and branch with **free** time and vector.

Every one-event signal is in $\mathcal M_1^{\rm TT}$, and so is every one-event signal affected by any per-arrival
timing error (statics, wrong medium) and any amplitude or polarization error. The **identifiability
margins** are

$$\delta_1(\xi)=\inf_{\boldsymbol\mu\in\mathcal M_1}\lVert\boldsymbol\mu_\xi-\boldsymbol\mu\rVert,\qquad
\delta_{\rm TT}(\xi)=\inf_{\boldsymbol\mu\in\mathcal M_1^{\rm TT}}\lVert\boldsymbol\mu_\xi-\boldsymbol\mu\rVert\ \le\ \delta_1(\xi).$$

Both are homogeneous of degree 1 in the amplitude: if the two tensors are multiplied by $\lambda$, the margin is
multiplied by $|\lambda|$. They are reported normalized by $\sqrt{E}$, with $E$ the network energy of one event
($\mathrm{SNR}=\sqrt E$ in units of noise standard deviation, the scale of Section 4 of the thesis).

## AB.A.1 The information limit

**Theorem AB.A.1 (Le Cam bound).** Let $\mathcal M$ be any set of signals and $\delta=\inf_{\boldsymbol\mu\in\mathcal M}\lVert\boldsymbol\mu_\xi-\boldsymbol\mu\rVert$.
For every test $\varphi$ of $H_0:\ \mathbf d=\boldsymbol\mu+\mathbf n,\ \boldsymbol\mu\in\mathcal M$ against
$H_1:\ \mathbf d=\boldsymbol\mu_\xi+\mathbf n$,

$$\sup_{\boldsymbol\mu\in\mathcal M}P_{\boldsymbol\mu}(\varphi=1)+P_{\xi}(\varphi=0)\ \ge\ 2\Phi(-\delta/2).$$

If the infimum is attained at $\boldsymbol\mu^\ast$, the Neyman-Pearson test between $\boldsymbol\mu^\ast$ and $\boldsymbol\mu_\xi$ with threshold at the
midpoint attains equality for that pair.

*Proof.* For fixed $\boldsymbol\mu\in\mathcal M$, $P_{\boldsymbol\mu}(\varphi=1)+P_\xi(\varphi=0)\ge1-\mathrm{TV}(\mathcal N(\boldsymbol\mu,I),\mathcal N(\boldsymbol\mu_\xi,I))$
(for any event $A$, $P(A^c)+Q(A)\ge1-\sup_A|P(A)-Q(A)|$). Projecting onto the direction
$\boldsymbol\mu_\xi-\boldsymbol\mu$, the total variation between two normals with identity covariance is that of two scalar normals
separated by $D=\lVert\boldsymbol\mu_\xi-\boldsymbol\mu\rVert$, that is $2\Phi(D/2)-1$. Thus the left-hand side is at least
$2\Phi(-D/2)$ for each $\boldsymbol\mu$, and the supremum is at least $\sup_{\boldsymbol\mu}2\Phi(-\lVert\boldsymbol\mu_\xi-\boldsymbol\mu\rVert/2)=2\Phi(-\delta/2)$
by continuity and monotonicity of $\Phi$. The equality for the midpoint test is the computation of its two errors,
$\Phi(-\delta/2)$ each. $\square$

**Corollary AB.A.1.** A total error $\le e$ requires $\delta\ge z_e:=2\Phi^{-1}(1-e/2)$; for $e=0.1$, $z_e=3.29$. With
$\mathcal M=\mathcal M_1^{\rm TT}$ the bound holds for any detector that must work even when the one-event
explanation has arbitrary per-arrival timing and amplitude errors.

What the theorem does not give: sufficiency for composite hypotheses. That a GLRT with unknown parameters attains the
limit up to constants is plausible and is not proved here.

## AB.A.2 Exact coincidence

**Proposition AB.A.2 (superposition).** If $\theta_1=\theta_2$, then $\boldsymbol\mu_\xi=G(\theta_1)(\mathbf m_1+\mathbf m_2)\in\mathcal M_1$ and
$\delta_1(\xi)=0$, whatever $\mathbf m_1$ and $\mathbf m_2$ are.

*Proof.* Linearity of $G(\theta)$ in the tensor. $\square$

The difference in mechanisms, however large, does not resolve two events at the same point and the same instant: the
sum of two tensors is a tensor. All the identifiability comes from the separation, and the difference in mechanisms
only changes **the order** with which the margin grows with it (Theorem AB.A.5).

## AB.A.3 The per-station relaxation is the timing-tolerant margin

**Theorem AB.A.3.** (a) $\delta_{\rm TT}(\xi)^2=\sum_{i=1}^{n_s}R_i$, with

$$R_i=\min_{t_P,t_S\in\mathbb R,\ \mathbf h_P,\mathbf h_S\in\mathbb R^3}\big\lVert\mathbf d_i-\mathbf h_P\bar s(\cdot-t_P)-\mathbf h_S\bar s(\cdot-t_S)\big\rVert^2
=\lVert\mathbf d_i\rVert^2-\max_{t_P,t_S}\ \frac{\lVert\mathbf b_P\rVert^2+\lVert\mathbf b_S\rVert^2-2\alpha\langle\mathbf b_P,\mathbf b_S\rangle}{1-\alpha^2},$$

where $\mathbf d_i$ is the part of $\boldsymbol\mu_\xi$ at station $i$, $\mathbf b_b=\langle\mathbf d_i,\bar s(\cdot-t_b)\rangle\in\mathbb R^3$ and
$\alpha=a(t_S-t_P)$ (at $|\alpha|=1$ the quotient is replaced by its limit, $\lVert\mathbf b_P\rVert^2$).
(b) $\delta_1\ge\delta_{\rm TT}$.
(c) If at each station the pulses of the P group (the two P of the data and the P of the model) are orthogonal to those of the
S group, $\delta_{\rm TT}$ does not change when the two events are shifted by the same statics $e_{ib}$ (common to both).

*Proof.* (a) The stations occupy orthogonal coordinates and $\mathcal M_1^{\rm TT}$ is a product of per-station
classes, so the squared distance is the sum of the per-station distances. For fixed times, the
least-squares fit of $(\mathbf h_P,\mathbf h_S)$ on two pulses with Gram matrix $\big(\begin{smallmatrix}1&\alpha\\\alpha&1\end{smallmatrix}\big)\otimes I_3$
leaves the residual $\lVert\mathbf d_i\rVert^2-\mathbf B^{\mathsf T}\big(\begin{smallmatrix}1&\alpha\\\alpha&1\end{smallmatrix}\big)^{-1}\mathbf B$, which is the expression
written. (b) $\mathcal M_1\subset\mathcal M_1^{\rm TT}$. (c) With orthogonal groups, the residual of each group is
invariant to a common translation of its data pulses, because the class allows the model pulse to be translated by the
same amount. $\square$

Part (c) is the useful content for a timing-tolerant method: the dominant model error (the station
statics: in the synthetic generator they are drawn with a standard deviation of 0.3 s, and 0.12 to 0.28 s remain after
the station corrections, thesis Section 6.2) shifts the two events equally and leaves intact the separations between
arrivals $\Delta_{ib}=t^2_{ib}-t^1_{ib}$, which is where identifiability comes from. Only the **event-specific**
errors (path differences between two nearby sources, of the order of $|\Delta\mathbf x|$ times the gradient of
the static) perturb the margin.

What the theorem does not give: at large separations the relaxation is loose, because the class can use its two free
pulses to fit the two S, the strongest ones. In the script, with $\Delta\tau=1$ s, $\delta_{\rm TT}=0.50$ and
$\delta_1=1.00$. For the overlap regime, which is the one that matters, the relaxation is almost exact (§AB.A.6).

## AB.A.4 One arrival: the exact formula

**Lemma AB.A.4.** For one arrival with two data pulses $\mathbf g_1\bar s+\mathbf g_2\bar s(\cdot-\Delta)$ and one model pulse,

$$R(\mathbf g_1,\mathbf g_2,\Delta)=\lVert\mathbf g_1\rVert^2+\lVert\mathbf g_2\rVert^2+2\langle\mathbf g_1,\mathbf g_2\rangle a(\Delta)-\max_{\varepsilon\in\mathbb R}\big\lVert\mathbf g_1a(\varepsilon)+\mathbf g_2a(\varepsilon-\Delta)\big\rVert^2 .$$

*Proof.* For fixed $\varepsilon$, the best $\mathbf h$ is $\langle\mathbf h_{\rm data},\bar s(\cdot-\varepsilon)\rangle$ because the pulse has
norm 1, and the residual is the energy minus $\lVert\mathbf h\rVert^2$. $\square$

## AB.A.5 Small separation: first order against second order

**Theorem AB.A.5.** Let $a\in C^6$ be even near 0 with the properties of the common notation, and $\mathbf G=\mathbf g_1+\mathbf g_2$.
As $\Delta\to0$ with $\mathbf g_1,\mathbf g_2$ fixed:

(i) if $\mathbf G\ne0$,
$$R=\omega^2\Delta^2\,\frac{\lVert\mathbf g_1\times\mathbf g_2\rVert^2}{\lVert\mathbf g_1+\mathbf g_2\rVert^2}+O(\Delta^4);$$

(ii) if $\mathbf g_1=\alpha\mathbf e$, $\mathbf g_2=\beta\mathbf e$ with $\lVert\mathbf e\rVert=1$ and $\alpha,\beta>0$ (parallel polarizations),
$$R=\frac{\psi-\omega^4}{4}\,\Delta^4\,\frac{\alpha^2\beta^2}{(\alpha+\beta)^2}+O(\Delta^6);$$

(iii) $\psi-\omega^4>0$ for every pulse in $L^2$.

*Proof.* We write $f(\varepsilon)=\lVert\mathbf g_1a(\varepsilon)+\mathbf g_2a(\varepsilon-\Delta)\rVert^2$. We will use three
global inequalities that come from $a(t)=\int|\hat s(\nu)|^2\cos(2\pi\nu t)\,d\nu$ and from
$1-y^2/2\le\cos y\le1-y^2/2+y^4/24$, $|\sin y|\le|y|$:

$$1-\tfrac{\omega^2t^2}{2}\le a(t)\le1-\tfrac{\omega^2t^2}{2}+\tfrac{\psi t^4}{24},\qquad |a'(t)|\le\omega^2|t| .\tag{A.1}$$

*Localization.* If $\mathbf G\ne0$: as $\Delta\to0$, $f\to\lVert\mathbf G\rVert^2a(\varepsilon)^2$ uniformly in $\varepsilon$
(because $a$ is Lipschitz), and $\sup_{|\varepsilon|\ge\eta}a(\varepsilon)^2<1$ for every $\eta>0$ (continuity,
$|a|<1$ away from 0 and $a\to0$). Since $f(0)\to\lVert\mathbf G\rVert^2$, the maximum is in $|\varepsilon|<\eta$ for small
$\Delta$. Inside, we write $\varepsilon=\Delta u$. With
$\mathbf g_1a(\Delta u)+\mathbf g_2a(\Delta(u-1))=\mathbf G\,a(\Delta u)+\mathbf g_2[a(\Delta(u-1))-a(\Delta u)]$, (A.1) gives for
$|\Delta u|\le\eta_0:=\min(\sqrt{6\omega^2/\psi},\,1/\omega)$ (where, by (A.1), $\tfrac12\le a(t)\le1-\omega^2t^2/4$)

$$f(\Delta u)\le\Big(\lVert\mathbf G\rVert\big(1-\tfrac{\omega^2\Delta^2u^2}{4}\big)+\lVert\mathbf g_2\rVert\,\omega^2\Delta^2(|u|+1)\Big)^2
=\lVert\mathbf G\rVert^2-\omega^2\Delta^2\Big(\tfrac{\lVert\mathbf G\rVert^2u^2}{2}-2\lVert\mathbf G\rVert\lVert\mathbf g_2\rVert(|u|+1)\Big)+O(\Delta^4u^4),$$

so there exists $U$ such that the maximum satisfies $|u|\le U$. *Expansion.* On $|u|\le U$, $a(x)=1-\omega^2x^2/2+O(x^4)$
gives uniformly

$$f(\Delta u)=\lVert\mathbf G\rVert^2-\omega^2\Delta^2q(u)+O(\Delta^4),\qquad q(u)=A u^2+B(u-1)^2,\quad A=\langle\mathbf g_1,\mathbf G\rangle,\ B=\langle\mathbf g_2,\mathbf G\rangle .$$

Since $A+B=\lVert\mathbf G\rVert^2>0$, $q$ is convex with minimum $AB/\lVert\mathbf G\rVert^2$ at $u^\ast=B/\lVert\mathbf G\rVert^2$, and
$\max f=\lVert\mathbf G\rVert^2-\omega^2\Delta^2AB/\lVert\mathbf G\rVert^2+O(\Delta^4)$. On the other hand
$\lVert\mathbf g_1\rVert^2+\lVert\mathbf g_2\rVert^2+2\langle\mathbf g_1,\mathbf g_2\rangle a(\Delta)=\lVert\mathbf G\rVert^2-\omega^2\Delta^2\langle\mathbf g_1,\mathbf g_2\rangle+O(\Delta^4)$.
Subtracting, $R=\omega^2\Delta^2\,(AB-\langle\mathbf g_1,\mathbf g_2\rangle\lVert\mathbf G\rVert^2)/\lVert\mathbf G\rVert^2+O(\Delta^4)$, and with
$c=\langle\mathbf g_1,\mathbf g_2\rangle$,
$AB-c\lVert\mathbf G\rVert^2=(\lVert\mathbf g_1\rVert^2+c)(\lVert\mathbf g_2\rVert^2+c)-c(\lVert\mathbf g_1\rVert^2+\lVert\mathbf g_2\rVert^2+2c)=\lVert\mathbf g_1\rVert^2\lVert\mathbf g_2\rVert^2-c^2=\lVert\mathbf g_1\times\mathbf g_2\rVert^2$
(Lagrange's identity). This proves (i).

(ii) The problem is scalar: $R=\alpha^2+\beta^2+2\alpha\beta a(\Delta)-\max_\varepsilon(\alpha a(\varepsilon)+\beta a(\varepsilon-\Delta))^2$, and
near 0 the inner function is positive, so it suffices to maximize it. With $S=\alpha+\beta$, $a(x)=1-\omega^2x^2/2+\psi x^4/24+O(x^6)$ and
the same localization, $\varepsilon=\Delta u$:

$$\alpha a(\Delta u)+\beta a(\Delta(u-1))=S-\tfrac{\omega^2\Delta^2}{2}\big(\alpha u^2+\beta(u-1)^2\big)+\tfrac{\psi\Delta^4}{24}\big(\alpha u^4+\beta(u-1)^4\big)+O(\Delta^6).$$

The order-2 term is maximized at $u^\ast=\beta/S$ with value $-\tfrac{\omega^2\Delta^2}{2}\tfrac{\alpha\beta}{S}$; moving $u$ by
$O(\Delta^2)$ around $u^\ast$ changes the value by $O(\Delta^6)$ because $u^\ast$ is stationary. At $u^\ast$,
$\alpha u^{\ast4}+\beta(u^\ast-1)^4=\alpha\beta(\alpha^3+\beta^3)/S^4$. Squaring,

$$\max(\cdot)^2=S^2-\omega^2\alpha\beta\Delta^2+\frac{\omega^4\alpha^2\beta^2\Delta^4}{4S^2}+\frac{\psi\Delta^4\alpha\beta(\alpha^3+\beta^3)}{12S^3}+O(\Delta^6),$$

while $\alpha^2+\beta^2+2\alpha\beta a(\Delta)=S^2-\omega^2\alpha\beta\Delta^2+\psi\alpha\beta\Delta^4/12+O(\Delta^6)$. The difference is
$\tfrac{\psi\alpha\beta\Delta^4}{12}\big(1-\tfrac{\alpha^3+\beta^3}{S^3}\big)-\tfrac{\omega^4\alpha^2\beta^2\Delta^4}{4S^2}$, and since
$S^3=\alpha^3+\beta^3+3\alpha\beta S$, the parenthesis equals $3\alpha\beta/S^2$. What remains is
$R=\tfrac{\alpha^2\beta^2\Delta^4}{4S^2}(\psi-\omega^4)+O(\Delta^6)$. The remainders are $O(\Delta^6)$ because $a$ is even.

(iii) $\langle\ddot{\bar s},\bar s\rangle=-\lVert\dot{\bar s}\rVert^2=-\omega^2$ (integration by parts), so
$\psi-\omega^4=\lVert\ddot{\bar s}\rVert^2\lVert\bar s\rVert^2-\langle\ddot{\bar s},\bar s\rangle^2\ge0$ by Cauchy-Schwarz, with
equality only if $\ddot{\bar s}\propto\bar s$, whose solutions are periodic and are not in $L^2$. $\square$

**Interpretation.** $\lVert\mathbf g_1\times\mathbf g_2\rVert=\lVert\mathbf g_1\rVert\lVert\mathbf g_2\rVert\sin\angle(\mathbf g_1,\mathbf g_2)$: what
resolves at first order is the **angle between the radiated amplitude vectors** in that arrival, that is, the difference
in mechanisms as seen by that station and that branch. Without it (parallel polarizations) the pulse of two coincident
components looks like a broadened pulse and only the curvature of the pulse, $\psi-\omega^4$, gives it away, at second
order. This is the angle "between subspaces at the stations where they overlap" requested by statement (a) of §6 of
the model design note, in quantitative form.

**Proposition AB.A.5' (large separation).** With $M=\max(\lVert\mathbf g_1\rVert,\lVert\mathbf g_2\rVert)$, $m=\min(\cdot)$ and
$\bar\alpha(\Delta)=\sup_{|u|\ge|\Delta|/2}|a(u)|$: $R\ge m^2(1-\bar\alpha^2)-4\bar\alpha Mm$.

*Proof.* For every $\varepsilon$, either $|\varepsilon|\ge|\Delta|/2$ or $|\varepsilon-\Delta|\ge|\Delta|/2$, so
$\lVert\mathbf g_1a(\varepsilon)+\mathbf g_2a(\varepsilon-\Delta)\rVert\le M+\bar\alpha m$; and the energy is at least
$M^2+m^2-2\bar\alpha Mm$. Subtracting the squares gives the bound. $\square$

With the 2 Hz Ricker, $\bar\alpha=1.8\times10^{-3}$ at 1.5 s and $1.2\times10^{-6}$ at 2 s: at those separations the weaker
arrival remains whole in the residual.

## AB.A.6 The network: geometry, mechanism and noise together

**Corollary AB.A.6.** Assume the P and S groups of each station are orthogonal (as in Theorem AB.A.3(c)) and
$\mathbf g^1_{ib}+\mathbf g^2_{ib}\ne0$ in every arrival. As the separation of the two events tends to 0,

$$\delta_{\rm TT}^2=\omega^2\sum_{i,b}\Delta_{ib}^2\,\kappa_{ib}+O(\max\Delta^4),\qquad \kappa_{ib}=\frac{\lVert\mathbf g^1_{ib}\times\mathbf g^2_{ib}\rVert^2}{\lVert\mathbf g^1_{ib}+\mathbf g^2_{ib}\rVert^2},$$

and with $\Delta_{ib}=\mathbf v_{ib}^{\mathsf T}(\theta_2-\theta_1)+O(|\theta_2-\theta_1|^2)$, $\mathbf v_{ib}=\nabla_\theta t_{ib}=(-\mathbf p_{ib},1)$ (the vector of AA.7(d)),

$$\delta_{\rm TT}^2\approx\omega^2\,(\theta_2-\theta_1)^{\mathsf T}F_{\rm pol}(\theta_2-\theta_1),\qquad F_{\rm pol}=\sum_{i,b}\kappa_{ib}\,\mathbf v_{ib}\mathbf v_{ib}^{\mathsf T}.$$

*Proof.* For small separations, at each station the optimal assignment of Theorem AB.A.3 puts one model pulse
in each group: any other leaves a whole group unexplained, with residual at least its energy, whereas
this assignment leaves $O(\Delta^2)$. With orthogonal groups $R_i=R_{iP}+R_{iS}$, and each term is the one of Theorem
AB.A.5(i). The linearization is the Taylor expansion of $t_{ib}$. $\square$

$F_{\rm pol}$ is the observability matrix of AA.8 with the energy weights $E_{ib}$ replaced by the
**polarization diversity** weights $\kappa_{ib}$. Three cases:

- **Co-located, different times.** $\Delta_{ib}=\Delta\tau$ for every arrival and
  $\delta_{\rm TT}^2=\omega^2\Delta\tau^2K_{\rm mec}+O(\Delta\tau^4)$, with $K_{\rm mec}=\sum_{ib}\kappa_{ib}$, the *mechanism diversity
  energy*, zero if and only if the radiated vectors are parallel in every arrival (for example, $\mathbf m_2\propto\mathbf m_1$).
  An upper bound of the same order comes from choosing the single event at $\theta_1$ with $\mathbf m=\mathbf m_1+\mathbf m_2$ and a common
  origin shift $\varepsilon$: $\delta_1^2\le\omega^2\Delta\tau^2\big[E_2-(\sum\langle\mathbf g^2,\mathbf G\rangle)^2/\sum\lVert\mathbf G\rVert^2\big]+O(\Delta\tau^3)$.
  Thus $\delta_1=\Theta(\Delta\tau)$ with explicit constants on both sides.
- **Co-located, same mechanism.** $K_{\rm mec}=0$ and Theorem AB.A.5(ii) governs: $\delta=\Theta(\Delta\tau^2)$. Since
  $\mathbf g^2_{ib}=r\,\mathbf g^1_{ib}$ with the same $r=\sqrt{E_2/E_1}$ in every arrival, the optimal shift $u^\ast$ is common, a
  single event realizes it, and $\delta_1=\delta_{\rm TT}(1+o(1))$.
- **Same mechanism, nearby positions.** $\kappa_{ib}=O(|\Delta\mathbf x|^2)$ because the radiation varies smoothly, and the two
  terms are $O(|\Delta\mathbf x|^4)$: $\delta=O(|\Delta\mathbf x|^2)$.

**Numerical verification** (`experiments/exploratory/THEORY-5/check_identifiability.py`, `experiments/exploratory/THEORY-5/check_identifiability_out.txt`).

*One arrival.* The ratio between the exact residual of Lemma AB.A.4 (maximized numerically) and the leading term of
Theorem AB.A.5 tends to 1: 0.907, 0.991, 0.999 at $\Delta=10,3,1$ ms with different polarizations, and 0.975, 0.997,
0.9998 at $\Delta=30,10,3$ ms with parallel polarizations. The bound of Proposition AB.A.5' holds in 200 random pairs
per separation (0.5 to 2 s).

*Network.* 8-station `spread` network, source at $(2,-3,9)$ km, two double couples of unit energy each. $\delta_1$ is
the exact distance, optimized over hypocentre, origin (Nelder-Mead from 5 starts) and tensor (least squares).

| Scenario | separation | exact $\delta_1/\sqrt E$ | $\delta_{\rm TT}/\sqrt E$ (Theorem AB.A.3) | Corollary AB.A.6 |
|:--|--:|--:|--:|--:|
| S1 co-located, different mechanisms | $\Delta\tau=3$ ms | 0.0235 | 0.0197 | 0.0197 |
| | 10 ms | 0.0784 | 0.0657 | 0.0657 |
| | 30 ms | 0.234 | 0.197 | 0.197 |
| | 100 ms | 0.729 | 0.427 | (0.657, out of range) |
| S2 co-located, same mechanism | $\Delta\tau=3$ ms | $2.809\times10^{-4}$ | $2.809\times10^{-4}$ | $2.809\times10^{-4}$ |
| | 30 ms | 0.0278 | 0.0278 | 0.0281 |
| | 100 ms | 0.278 | 0.276 | (0.312) |
| S3 same mechanism and origin | $\Delta x=0.1$ km | 0.0127 | 0.0098 | 0.0098 |
| | 0.3 km | 0.108 | 0.080 | 0.088 |
| | 1 km | 0.678 | 0.496 | (0.93) |

The log-log slopes at small separation are **1.00** (S1), **2.00** (S2) and **1.97** (S3), those of Theorem AB.A.5. The
timing-tolerant relaxation is 84 % of the exact margin in S1, 100 % in S2 and 77 % in S3: giving up the
predicted per-arrival times costs little near overlap. In all rows $\delta_{\rm TT}\le\delta_1$, as
required by Theorem AB.A.3(b) (the first version of the script violated this inequality because of a badly initialized
local optimizer; the stored version uses a vectorized global search and satisfies it in all rows).

## AB.A.7 Minimum resolvable separation

**Corollary AB.A.7.** With Theorem AB.A.1 and Corollary AB.A.6, and amplitudes scaled to SNR $\sqrt E$ per event, every
detector with total error $\le e$ needs, at small separation,

$$\Delta\tau\ \gtrsim\ \frac{z_e}{\omega\sqrt{K_{\rm mec}}}\ \propto\ \mathrm{SNR}^{-1}\quad\text{(different mechanisms)},\qquad
\Delta\tau\ \gtrsim\ \Big(\frac{2z_e}{\sqrt{\psi-\omega^4}\,\sqrt{\sum_{ib}\big(\tfrac{|\mathbf g^1_{ib}||\mathbf g^2_{ib}|}{|\mathbf g^1_{ib}|+|\mathbf g^2_{ib}|}\big)^2}}\Big)^{1/2}\ \propto\ \mathrm{SNR}^{-1/2}\quad\text{(same mechanism)}.$$

These are **necessary** conditions (information limit with the exact model). The second follows from
$\delta_{\rm TT}^2=\tfrac{\psi-\omega^4}{4}\Delta\tau^4\sum_{ib}(\cdots)^2$ (Theorem AB.A.5(ii) summed over arrivals)
and $\delta\ge z_e$; the version first circulated had $4z_e$ in place of $2z_e$, a factor $\sqrt2$ too large in
$\Delta\tau$ (corrected on integration; the scaling and the table below, which comes from the exact curve, are
unaffected). Interpolating the exact curve $\delta_1$ of the
script, for total error 0.1:

| Scenario | SNR 12.4 (2M threshold) | SNR 30 | SNR 100 |
|:--|--:|--:|--:|
| S1 co-located, different mechanisms | 0.034 s | 0.014 s | 0.004 s |
| S2 co-located, same mechanism | 0.098 s | 0.061 s | 0.033 s |
| S3 same mechanism and origin, separated in plan view | 0.54 km | 0.30 km | 0.16 km |

Reading for the atlas (roadmap, phase 5, primary line): the "separation between arrivals" axis of the atlas grid
(< 0.5, 0.5–1, 1–3, > 3 s, §7 of the model design note) lies entirely above the information limit for
SNR $\ge12$; the < 0.5 s cell should be split into < 0.1 s (where identifiability starts to dominate, above all with
similar mechanisms) and 0.1–0.5 s. The covariate that the atlas must record, besides the separation, is
$K_{\rm mec}/E$ (computable for each injected pair), because it changes the order of the scaling law.

**Conjecture AB.A.8 (Fisher information of the pair).** The Fisher matrix of the 20 parameters of the pair becomes singular
as $\Delta\tau\to0$, with its smallest eigenvalue of order $\Delta\tau^2$ if $K_{\rm mec}>0$ and of order $\Delta\tau^4$ if
$K_{\rm mec}=0$, and the Cramér-Rao bound of $\Delta\tau$ scales accordingly. It is the estimation version of AB.A.5; it is not
proved here.

---

# Part B. Greedy recovery with a dictionary perturbed by timing error

## AB.B.0 The algorithm and the perturbation

The matched-filter record is $\mathbf y_i(t)=\sum_k\sum_b\mathbf g^k_{ib}\,a(t-t^k_{ib})+\boldsymbol\zeta_i(t)$, with
$\boldsymbol\zeta$ the filtered noise and the true times $t^k_{ib}=t_{ib}(\theta_k)+e_{ib}+e^k_{ib}$: common statics
$e_{ib}$ and event-specific errors $e^k_{ib}$. The dictionary is the grid nodes $j$, with predicted times
$t_{ib}(\theta_j)$, windows $I_{ib}(j)=\lfloor t_{ib}(\theta_j)/\Delta\rfloor+\{-W,\dots,W\}$ and projectors $\Pi_j$ of rank 6. The
deployed detector (MT-FW without cap, equal to `fw_record` with infinite cap; with cap, see Remark AB.B.7):

1. reads in each block the highest-energy vector of the window, $\mathbf v_{j,ib}=\mathbf y_i(\hat t_{ib})$, $\hat t_{ib}=\arg\max_{t\in I_{ib}(j)}\lVert\mathbf y_i(t)\rVert$ (D1);
2. scores $T_j=\lVert\Pi_j\mathbf v_j\rVert$ and chooses the maximum among the admissible nodes;
3. fits the tensor by least squares to $\mathbf v_j$ and subtracts $\mathbf p_{ib}\,a(t-\hat t_{ib})$, with $\mathbf p=\Pi_j\mathbf v_j$, at the
   **read** times, which are on the sample grid;
4. excludes the nodes within 5 km of the accepted one and repeats.

"Perturbed dictionary" means three things at once here: the true arrivals are shifted with respect to the
predicted ones (statics and event-specific error), they are off the sample grid (offset $\le\Delta/2$ with respect to the
nearest sample), and, if the source is not at a node, the subspace of the nearest node does not exactly contain
the signal (mismatch $\eta_k=\lVert(I-\Pi_{j_k})\hat{\mathbf g}_k\rVert$).

For source $k$ we call **representative node** $j_k$ one whose windows contain its arrivals, and **hit
set** $H_k$ the nodes within 5 km and 1 s (the repository's LOC rule). The rest of the admissible nodes form the
**wrong** set $\mathcal W$.

## AB.B.1 When D1 reads the correct arrival

**Lemma AB.B.1 (dominance of the peak reading).** In one block, let $\mathbf y(t)=\mathbf g\,a(t-t_0)+\boldsymbol\zeta(t)$ in the
window $I$, let $t^\circ\in I$ be the sample nearest to $t_0$, with offset $\delta^\circ=|t^\circ-t_0|\le\Delta/2$, and
$Z=\max_{t\in I}\lVert\boldsymbol\zeta(t)\rVert$. If

$$\lambda:=a(\delta^\circ)-\frac{2Z}{\lVert\mathbf g\rVert}\ >\ a^-=0.618,$$

then the sample read satisfies $a(\hat t-t_0)\ge\lambda$ and the vector read is $\mathbf v=a(\hat t-t_0)\,\mathbf g+\boldsymbol\zeta(\hat t)$
with $a(\hat t-t_0)\in[\lambda,1]$. In particular $\hat t$ is in the main lobe, within
$a_+^{-1}(\lambda)$ of $t_0$ ($a_+^{-1}$ the inverse of $a$ on $[0,0.118\ \mathrm s]$).

*Proof.* $\lVert\mathbf y(\hat t)\rVert\ge\lVert\mathbf y(t^\circ)\rVert\ge\lVert\mathbf g\rVert a(\delta^\circ)-Z$ by the choice
of the maximum, and $\lVert\mathbf y(\hat t)\rVert\le\lVert\mathbf g\rVert|a(\hat t-t_0)|+Z$. Thus $|a(\hat t-t_0)|\ge\lambda>a^-$; if
$a(\hat t-t_0)$ were negative, its absolute value would be at most $a^-$, a contradiction. $\square$

At 16 Hz, $a(\Delta/2)=0.906$, so the condition is $Z<0.144\,\lVert\mathbf g\rVert$ at the worst offset. The interference
$Z$ includes the noise and **the arrivals of another event within the same window**: if a second arrival of amplitude
$\lVert\mathbf g'\rVert$ falls at $\Delta'$ from the first, it contributes $\lVert\mathbf g'\rVert|a(\Delta')|$ to $Z$. The negative lobe of
$-0.618$ at 0.216 s, inside the tolerance of $\pm0.25$ s, is what makes D1 jump arrivals when two are
close: it is the quantitative version of the "wrong peak" of R1 (36 % of P and 20 % of S at more than 0.5 s) that the
model design note identifies as the root cause.

## AB.B.2 Lower bound at the representative node

**Lemma AB.B.2.** Let $D$ be the set of blocks of node $j_k$ that satisfy Lemma AB.B.1 (with their $\lambda_{ib}$ and their $Z_{ib}$),
$Y_{ib}=\max_{t\in I_{ib}}\lVert\mathbf y_i(t)\rVert$ and $\lVert\mathbf Y\rVert=(\sum_{ib}Y^2_{ib})^{1/2}$. Then

$$T_{j_k}\ \ge\ \mathrm{LB}_k:=\frac1{\lVert\mathbf g_k\rVert}\Big[\sum_{ib\in D}\big(\lambda_{ib}\lVert\mathbf g_{ib}\rVert^2-\lVert\mathbf g_{ib}\rVert Z_{ib}\big)-\sum_{ib\notin D}\lVert\mathbf g_{ib}\rVert Y_{ib}\Big]-\eta_k\lVert\mathbf Y\rVert .$$

*Proof.* For $\mathbf u=\Pi_{j_k}\hat{\mathbf g}_k$, $\lVert\mathbf u\rVert\le1$ and $T_{j_k}=\lVert\Pi_{j_k}\mathbf v\rVert\ge\langle\mathbf u,\mathbf v\rangle$
when the latter is positive. $\langle\mathbf u,\mathbf v\rangle=\langle\hat{\mathbf g},\mathbf v\rangle-\langle(I-\Pi)\hat{\mathbf g},\mathbf v\rangle\ge\langle\hat{\mathbf g},\mathbf v\rangle-\eta_k\lVert\mathbf Y\rVert$.
In the blocks of $D$, by Lemma AB.B.1,
$\langle\hat{\mathbf g}_{ib},\mathbf v_{ib}\rangle\ge(\lambda_{ib}\lVert\mathbf g_{ib}\rVert^2-\lVert\mathbf g_{ib}\rVert Z_{ib})/\lVert\mathbf g\rVert$; outside,
$|\langle\hat{\mathbf g}_{ib},\mathbf v_{ib}\rangle|\le\lVert\mathbf g_{ib}\rVert Y_{ib}/\lVert\mathbf g\rVert$. $\square$

A block whose true arrival falls outside the window (timing error larger than the tolerance) is outside $D$ and
**subtracts**: insufficient tolerance does not only lose energy, it turns it into noise working against the detection.

## AB.B.3 Leakage towards a wrong node

**Lemma AB.B.3.** For a node $j$ and a source $k$, let $\alpha^{jk}_{ib,b'}=\max_{t\in I_{ib}(j)}|a(t-t^k_{ib'})|$ (how much of
arrival $b'$ of source $k$ at station $i$ the window $(i,b)$ of node $j$ can read) and let
$\mathbf W_r\in\mathbb R^{3n_b}$, $r=(ib,b')$, be the vectors with $\alpha^{jk}_{ib,b'}\mathbf g^k_{ib'}$ in block $ib$ and zero outside. For
every choice of reading samples, the signal part of the reading satisfies

$$\lVert\Pi_j\mathbf v^{\rm signal}_j\rVert\le\sum_k\Lambda_{jk},\qquad \Lambda_{jk}=\max_{\boldsymbol\sigma\in\{\pm1\}^r}\Big\lVert\Pi_j\sum_r\sigma_r\mathbf W_r\Big\rVert
\ \le\ \Big(\sum_{ib}\big(\textstyle\sum_{b'}\alpha^{jk}_{ib,b'}\lVert\mathbf g^k_{ib'}\rVert\big)^2\Big)^{1/2},$$

and the noise part satisfies $\lVert\Pi_j\mathbf v^{\rm noise}_j\rVert\le N_j:=\max_{\mathbf s}\lVert\Pi_j\boldsymbol\zeta_{\mathbf s}\rVert$ (maximum over
all choices of samples). Under Gaussian noise, with $N$ nodes,
$P\big(\max_jN_j\ge\sqrt6+\sqrt{2n_b\log(2W+1)}+\sqrt{2\log(N/p)}\big)\le p$ (Proposition AB.C.1).

*Proof.* $\lVert\Pi_j\mathbf x\rVert=\max_{\mathbf u\in\mathrm{ran}\,\Pi_j,\lVert\mathbf u\rVert\le1}\langle\mathbf u,\mathbf x\rangle$. The signal read in
block $ib$ at sample $t$ is $\sum_k\sum_{b'}\mathbf g^k_{ib'}a(t-t^k_{ib'})$, so
$\langle\mathbf u,\mathbf v^{\rm signal}\rangle\le\sum_k\sum_r\alpha_r|\langle\mathbf u_{ib},\mathbf g^k_{ib'}\rangle|=\sum_k\max_{\boldsymbol\sigma}\langle\mathbf u,\sum_r\sigma_r\mathbf W^k_r\rangle$;
the maximum over $\mathbf u$ of a sum is at most the sum of the maxima. The Cauchy-Schwarz form bounds each
$|\langle\mathbf u_{ib},\cdot\rangle|$ by $\lVert\mathbf u_{ib}\rVert\lVert\cdot\rVert$. The noise bound is Proposition AB.C.1 applied to each
node (the maximum over $\mathbf s$ of 1-Lipschitz functions is 1-Lipschitz) and a union over nodes. $\square$

$\Lambda_{jk}$ is the **timing-tolerant coherence** between source $k$ and node $j$: it is computed exactly
by enumerating $2^{r-1}$ signs (in practice $r\le16$; the leftover terms are bounded by Cauchy-Schwarz). Unlike
the cell bound AA.6, it sees the geometry of the subspace.

## AB.B.4 Correct selection and recovery

**Theorem AB.B.4 (one step).** Let $\mathcal R$ be the sources not yet found and assume that their representative nodes are
admissible. If

$$\max_{k\in\mathcal R}\mathrm{LB}_k\ >\ \max_{j\in\mathcal W}\Big[\sum_{k\in\mathcal R}\Lambda_{jk}+\Lambda^{\rm rest}_j\Big],$$

with $\Lambda^{\rm rest}_j$ any bound on $\lVert\Pi_j(\cdot)\rVert$ for the reading of the rest of the record (noise plus residuals
of sources already found), then the chosen hypothesis is in $\bigcup_{k\in\mathcal R}H_k$.

*Proof.* $T_{\rm chosen}\ge T_{j_k}\ge\mathrm{LB}_k$ (Lemma AB.B.2) and $T_j\le\sum_k\Lambda_{jk}+\Lambda^{\rm rest}_j$ for every
$j\in\mathcal W$ (Lemma AB.B.3 and the triangle inequality). $\square$

**Theorem AB.B.5 ($K$ steps).** If at each step $r=1,\dots,K$ the condition of Theorem AB.B.4 holds with $\mathcal R$ the sources
not found, $\mathcal W_r$ the admissible nodes outside $\bigcup_{k\in\mathcal R}H_k$ (the hit nodes of sources already
found that remain outside the exclusion count as wrong, so that there are no duplicates) and $\Lambda^{\rm rest}$
including the deflation residuals (Proposition AB.B.6), then the $K$ selections hit $K$ distinct
sources.

*Proof.* Induction on $r$ with Theorem AB.B.4. $\square$

**Proposition AB.B.6 (what the deployed deflation leaves).** Isolated source, node on the grid ($\eta=0$), no noise,
reading at the nearest sample with offset $\delta_{ib}\le\Delta/2$, so that $\mathbf v_{ib}=a(\delta_{ib})\mathbf g_{ib}$. After
subtracting $\mathbf p_{ib}a(t-\hat t_{ib})$ with $\mathbf p=\Pi_j\mathbf v$, the residual at station $i$ and branch $b$ is

$$\mathbf L_{ib}(t)=\mathbf g_{ib}\big[a(t-t_0)-a(\delta_{ib})\,a(t-t_0-\delta_{ib})\big]-\big(\mathbf p-\boldsymbol\lambda\circ\mathbf g\big)_{ib}\,a(t-\hat t_{ib}),$$

with $\boldsymbol\lambda\circ\mathbf g$ the vector of blocks $a(\delta_{ib})\mathbf g_{ib}$, and

$$\sup_t\lVert\mathbf L_{ib}(t)\rVert\ \le\ \lVert\mathbf g_{ib}\rVert\,c(\delta_{ib})+\lVert(I-\Pi_j)(\boldsymbol\lambda\circ\mathbf g)\rVert,\qquad c(\delta)=\sup_x|a(x)-a(\delta)a(x-\delta)|,$$

$$\lVert(I-\Pi_j)(\boldsymbol\lambda\circ\mathbf g)\rVert\le\min_{\bar\lambda}\Big(\sum_{ib}(\lambda_{ib}-\bar\lambda)^2\lVert\mathbf g_{ib}\rVert^2\Big)^{1/2}\le\tfrac{1-a(\Delta/2)}2\lVert\mathbf g\rVert .$$

*Proof.* The decomposition is algebra with $\hat t_{ib}=t_0+\delta_{ib}$. The first term is bounded by the definition of
$c$. For the second, $\mathbf p-\boldsymbol\lambda\circ\mathbf g=-(I-\Pi_j)(\boldsymbol\lambda\circ\mathbf g)$, and since $\bar\lambda\mathbf g\in\mathrm{ran}\,\Pi_j$,
$(I-\Pi_j)(\boldsymbol\lambda\circ\mathbf g)=(I-\Pi_j)((\boldsymbol\lambda-\bar\lambda)\circ\mathbf g)$, with norm at most that of its argument; with
$\lambda_{ib}\in[a(\Delta/2),1]$ and $\bar\lambda$ the midpoint the last bound follows. $\square$

The deployed deflation subtracts the waveform at the sample read, not at the true time. That is why, without noise and with everything
exact, it leaves in the arrivals of a found source a fixed fraction of its amplitude, which does not decrease with the SNR. For a
strong event that residual is far above the noise and is interference ($Z$ in Lemma AB.B.1) for any other
overlapping event. It is a "ghost" mechanism different from the coda one pointed out by the model design note (H4),
and it is corrected without changing the dictionary: subtracting at the refined sub-sample time (interpolation of the matched-filter
peak) and with the amplitude refitted at that time cancels the first term.

**Remark AB.B.7 (with cap and with ACAP).** With a radial map $\phi$, Lemmas AB.B.2 and AB.B.3 hold with $\phi(\mathbf v)$ in place
of $\mathbf v$; the cap applied to blocks with $\lVert\mathbf v_{ib}\rVert>c$ changes $\lambda_{ib}$ to $\lambda_{ib}\,c/\lVert\mathbf v_{ib}\rVert$ and the
residual of AB.B.6 adds the term of AA.4(g). ACAP does not alter the blocks when max/median $\le c_1$ (Proposition AB.C.4(c)).

**Remark AB.B.8 (the trade-off of the tolerance $W$).** $W\Delta$ must cover the timing error (plus the offset of the
node) for the blocks of the source to be in $D$; each extra sample increases the noise bound as
$\sqrt{2n_b\log(2W+1)}$ (with $n_b=16$: 14.9, 16.1 and 17.3 for $W=2,4,8$, $p=0.05$ and the grid of the script) and, above all,
the leakage $\Lambda_{jk}$, because more wrong nodes reach the true arrivals and the negative lobe of the
autocorrelation with their windows. The optimum is the smallest tolerance that covers the timing error, which supports
narrowing the tolerance where the station terms are reliable (the Gaussian process of §4.2 of the model design
note, which gives a local tolerance).

## AB.B.5 Numerical verification

`experiments/exploratory/THEORY-5/check_greedy.py`, `experiments/exploratory/THEORY-5/check_greedy_out.txt`: the deployed detector (without cap) on the
repository grid (3 km, 15 × 15 × 3 nodes, origins every 0.1 s, 61,425 hypotheses), two double-couple sources at nodes,
separated by at least 9 km in plan view, with statics $e_{ib}\sim U[-\varepsilon,\varepsilon]$ common to both, arrivals off the
sample grid and filtered Gaussian noise; 12 replicates per cell.

Results (`check_greedy_out.txt`, and `diag_greedy.py` / `diag_greedy_out.txt` for the slack):

- **The lemmas always hold.** The lower bound of Lemma AB.B.2 never exceeded the actual score of the true node (0
  of 432 cases), the leakage of Lemma AB.B.3 never fell below the noise-free score of a wrong node (0 of 191), and
  no certificate accompanied a wrong selection.
- **The certificate never activated**, whereas the algorithm was correct at the first step in 100 % of the 216
  replicates (SNR 30, 60 and 120 per event; $\Delta\tau=0$, 1 and 3 s; $\varepsilon=0$ and 0.15 s) and found both sources in all
  cells but one (0.92 with SNR 60, $\Delta\tau=3$ s, $\varepsilon=0.15$ s). Theorem AB.B.4 is correct but, in this geometry,
  vacuous.
- **Why** (diagnosis with SNR 120, $\Delta\tau=3$ s, $\varepsilon=0.15$ s, 8 replicates; medians relative to the score of the true
  node): the lower bound stays at 0.59; the best **actual wrong** node reaches 0.85 (at 3 to 6 km and 1 s from the
  true one in most cases, and in one replicate at 1.03); and the leakage-plus-noise bound rises to 1.11. The actual
  margin between the correct hypothesis and the best wrong one is of the order of 15 %, and the relaxations by signs and by
  window maxima lose more than that. The useful finding is not the certificate but the figure: **with statics of
  0.15 s the score surface is almost flat just outside the hit rule**, so that small changes
  (another LOC rule, another grid, another tolerance) move hits and misses. This agrees with the sensitivity of the
  headline figures to the hit rule that the audit noted (LOC against LOC+T).
- **The deflation residual is the one predicted by Proposition AB.B.6.** In the blocks with $\lVert\mathbf g_{ib}\rVert\ge0.1$ of the
  maximum, the maximum residual within $\pm0.5$ s of the arrival is a median fraction of **0.16**, 90th percentile of 0.32 and maximum of
  0.36 of the amplitude of the arrival, against $c(\delta)=0.077,0.154,0.232,0.310,0.387$ for offsets of 6, 12, 19, 25 and 31 ms. An
  event of 100 noise standard deviations at a station leaves there, without noise and with the exact model, a residual of 16 to
  36 standard deviations: it is interference for any other overlapping event and a source of duplicate detections.

To make the theorem useful in the overlap regime two changes are needed, which remain as future work: taking as the
lower bound the maximum over all of $H_k$, not only the representative node, and replacing the relaxation by signs with the
exact leakage with joint reading selection (§4.3 of the model design note), which is the one that sees the subspace
difference between neighbouring nodes.

---

# Part C. The adaptive robust statistic

## AB.C.0 Setting

That of AA.9: at hypothesis $h$, $n_b=2n_s$ blocks $\mathbf v_k\in\mathbb R^3$ are read, a map $\phi$ is applied block by block and
the score is $\sigma_\phi(h)=\lVert Q_h^{\mathsf T}\phi(\mathbf v)\rVert$, with $Q_h\in\mathbb R^{3n_b\times6}$ with orthonormal columns and blocks
$Q_k\in\mathbb R^{3\times6}$. Three properties are required: bounded influence, non-saturation and a false-alarm tail tighter
than that of AA.9.5.

## AB.C.1 Correction of the bound of AA.9.5(b)

**Proposition AB.C.1.** Under the assumptions of AA.9.5 (filtered white Gaussian noise with unit variance, independent components
and blocks, peak reading among $L=2W+1$ samples) and for $\phi$ the identity, the cap or the
logarithmic map,

$$P\Big(\sigma_\phi\ge\sqrt6+\sqrt{2n_b\log L}+\sqrt{2\log(1/p)}\Big)\le p .$$

That is, the constant $\kappa=\max_f|\hat k(f)|/\lVert k\rVert_2=2.146$ of AA.9.5(b) can be replaced by 1 in both
terms.

*Proof.* For a fixed choice $\mathbf s$ of one sample per block, the linear map $\mathbf z\mapsto\mathbf x_{\mathbf s}$ from the
white noise to the $3n_b$ samples read has one row per (block, component) pair; each row is the shifted normalized
kernel, of norm 1, and distinct rows have disjoint supports in $\mathbf z$ (distinct components,
distinct stations, or P and S of one station separated by more than the support of the kernel, which is the independence
assumption of AA.9.5). It is a partial isometry, of norm 1. With the 1-Lipschitz radial map and the contraction
of Lemma AA.9.1, $f_{\mathbf s}=\lVert Q_h^{\mathsf T}\phi(\mathbf x_{\mathbf s})\rVert$ is 1-Lipschitz in $\mathbf z$, and so is $F=\max_{\mathbf s}f_{\mathbf s}\ge\sigma_\phi$.
Gaussian concentration gives $P(F\ge\mathbb EF+t)\le e^{-t^2/2}$. Each $f_{\mathbf s}-\mathbb Ef_{\mathbf s}$ is
sub-Gaussian with parameter 1 (again by concentration) and $\mathbb Ef_{\mathbf s}\le(\mathbb Ef_{\mathbf s}^2)^{1/2}\le\sqrt6$ as in AA.9.5,
so the maximum of the $L^{n_b}$ satisfies $\mathbb EF\le\sqrt6+\sqrt{2n_b\log L}$. $\square$

The constant $\kappa$ is the norm of the full convolution, the correct one for a function of the **whole** record; the
reading uses one sample per block and component, and for that the constant is 1. The bound of AA.9.5(b) is not false,
only loose; the thesis (Section 4) and AA.9.5 now cite the new one, and the "factor of four" became a "factor of two".

| $p$ | AA.9.5(b) ($\kappa$) | Proposition AB.C.1 | empirical quantile of $\sigma$ (identity, cap, log, ACAP) |
|--:|--:|--:|:--|
| $10^{-2}$ | 25.80 | 13.33 | 5.82, 5.82, 5.52, 5.82 |
| $10^{-3}$ | 27.26 | 14.01 | 6.62, 6.62, 6.26, 6.62 |
| $10^{-4}$ | 28.50 | **14.59** | 7.12, 7.12, 6.76, 7.12 |

($2\times10^5$ hypotheses, same configuration as `code/theory_tail_check.py`: $n_b=14$, $W=4$; `experiments/exploratory/THEORY-5/check_adaptive.py`.
The quantiles reproduce those of `results/theory/tail_check_out.txt` with another seed: 7.18 and 6.78 there.) What remains
of the slack is the union over the $L^{n_b}$ possible readings.

## AB.C.2 A conditional tail that does not depend on the law of the amplitudes

**Assumption (I).** Under $H_0$, the read blocks $\mathbf v_1,\dots,\mathbf v_{n_b}$ are independent and the law of each one is
invariant under rotations of $\mathbb R^3$. The laws of the norms are arbitrary: heavy tails, different per block, non-stationary,
with transients of any amplitude provided that their direction is isotropic.

Under Gaussian noise the assumption holds with the peak reading: the $L\times3$ sample matrix of a block has a law
invariant under $X\mapsto XR$ (the three components are independent copies of the same process) and the energy $\arg\max$
is invariant under rotations, so $\mathbf v\mapsto R\mathbf v$ preserves the law.

**Norm-adaptive maps.** $\phi(\mathbf v)_k=\rho_k\,\mathbf v_k/\lVert\mathbf v_k\rVert$ with $\rho_k=r_k(\lVert\mathbf v_1\rVert,\dots,\lVert\mathbf v_{n_b}\rVert)\ge0$:
this includes the identity, the cap, the logarithmic map and ACAP (§AB.C.4).

**Theorem AB.C.2.** Under (I), for every norm-adaptive map, with $M=\sum_k\rho_k^2Q_k^{\mathsf T}Q_k\in\mathbb R^{6\times6}$ and every $t>0$,

$$P\Big(\sigma_\phi^2>\tfrac13\big[\mathrm{tr}\,M+2\sqrt{\mathrm{tr}(M^2)\,t}+2\lVert M\rVert t\big]\ \Big|\ \lVert\mathbf v_1\rVert,\dots,\lVert\mathbf v_{n_b}\rVert\Big)\le e^{-t},
\qquad \mathbb E\big[\sigma_\phi^2\,\big|\,\text{norms}\big]=\tfrac13\mathrm{tr}\,M .$$

*Proof.* (1) A rotation-invariant law on $\mathbb R^3$ has a direction uniform on $S^2$ and independent of the
norm (polar decomposition of the measure); with independent blocks, the directions $\boldsymbol\theta_k$ are i.i.d.
uniform and independent of all the norms, and therefore of the $\rho_k$. (2) $Q_h^{\mathsf T}\phi(\mathbf v)=A\mathbf x$ with
$\mathbf x=(\boldsymbol\theta_1,\dots,\boldsymbol\theta_{n_b})\in\mathbb R^{3n_b}$ and $A=Q_h^{\mathsf T}\mathrm{diag}(\rho_kI_3)$, fixed given the norms. (3) $\mathbf x$ is
sub-Gaussian with parameter $1/3$: by Archimedes' theorem, $\langle\boldsymbol\alpha,\boldsymbol\theta\rangle$ is uniform on
$[-\lVert\boldsymbol\alpha\rVert,\lVert\boldsymbol\alpha\rVert]$, so $\mathbb Ee^{\langle\boldsymbol\alpha,\boldsymbol\theta\rangle}=\sinh\lVert\boldsymbol\alpha\rVert/\lVert\boldsymbol\alpha\rVert$, and
$\sinh y/y=\sum_jy^{2j}/(2j+1)!\le\sum_j(y^2/6)^j/j!=e^{y^2/6}$ because $(2j+1)!\ge6^jj!$ (induction:
$(2j+3)(2j+2)\ge6(j+1)$); with independence, $\mathbb Ee^{\langle\boldsymbol\alpha,\mathbf x\rangle}\le e^{\lVert\boldsymbol\alpha\rVert^2/6}$. (4) Theorem 1
of Hsu, Kakade and Zhang (2012) for quadratic forms of centred sub-Gaussian vectors gives
$P(\lVert A\mathbf x\rVert^2>\sigma^2[\mathrm{tr}\,\Sigma+2\sqrt{\mathrm{tr}(\Sigma^2)t}+2\lVert\Sigma\rVert t])\le e^{-t}$ with $\Sigma=A^{\mathsf T}A$ and
$\sigma^2=1/3$; the non-zero spectrum of $A^{\mathsf T}A$ is that of $AA^{\mathsf T}=\sum_k\rho_k^2Q_k^{\mathsf T}Q_k=M$. The mean follows from
$\mathbb E\boldsymbol\theta\boldsymbol\theta^{\mathsf T}=I_3/3$. $\square$

**Corollary AB.C.2 (conditional threshold and p-value).** The threshold
$\tau_p(h)=\big(\tfrac13[\mathrm{tr}M+2\sqrt{\mathrm{tr}(M^2)\log(1/p)}+2\lVert M\rVert\log(1/p)]\big)^{1/2}$ has false-alarm
probability $\le p$ per hypothesis, conditional and therefore unconditional, **for any law of the amplitudes**; and
$\hat p(h)=e^{-t^\ast}$, with $\sqrt{t^\ast}=\big(-\sqrt{\mathrm{tr}M^2}+\sqrt{\mathrm{tr}M^2-2\lVert M\rVert(\mathrm{tr}M-3\sigma_\phi^2)}\big)/(2\lVert M\rVert)$
when $3\sigma^2_\phi>\mathrm{tr}M$ (and $\hat p=1$ otherwise), is a valid p-value in closed form. Over $N$ hypotheses, $p/N$.

$M$ is computed with the read norms, which the detector already has: the threshold adapts by itself to the background level at each
station, each hour and each hypothesis. $\mathrm{tr}M=\sum_k\rho_k^2\lVert Q_k\rVert_F^2$ is the **incoherent** energy (the
expectation of the projection if the directions were random) and $\sigma_\phi^2$ the **coherent** one; the test compares the two,
like a semblance with the subspace of the tensor as the coherence model.

**Verification** (`experiments/exploratory/THEORY-5/check_adaptive_out.txt`):

| | identity | cap | log | ACAP |
|:--|--:|--:|--:|--:|
| Gaussian: median conditional threshold at $p=10^{-4}$ | 9.75 | 9.75 | 9.24 | 9.75 |
| Gaussian: realized rate at $p=10^{-2},10^{-3},10^{-4}$ ($2\times10^5$ hypotheses) | 0 | 0 | 0 | 0 |
| Heavy tails + isotropic transients: realized rate at $p=10^{-3}$ ($5\times10^4$) | 0 | 0 | 0 | 0 |
| same, with the fixed Gaussian threshold of $10^{-3}$ (readings divided by the true scale) | 0.96 | 0.54 | 0.81 | 0.89 |

The conditional threshold lies within a factor 1.37 of the unconditional empirical quantile (7.12) in units of $\sigma$, against
2.0 for Proposition AB.C.1 and 4.0 for AA.9.5(b); against the exact conditional quantile (Monte Carlo over directions, 400
cases, $p=10^{-3}$) the median ratio is 1.40 (1.33 to 1.55). The heavy-tailed noise is an isotropic Student $t$
with $\nu=2$ per sample, a lognormal scale ($\sigma=1$) per block and transients in 15 % of the blocks with
log-uniform amplitude between $10^2$ and $10^4$. There any fixed threshold fails (rates of 0.54 to 0.96 where
$10^{-3}$ was required) and the conditional one recorded no false alarm.

Limits of Theorem AB.C.2, stated in advance:

- **Polarization.** Real noise is not isotropic: the horizontal and vertical power differ, and there are polarized
  transients. For the identity, if $\mathbf v_k=\Lambda_k^{1/2}\mathbf w_k$ with $\mathbf w_k$ isotropic and $\Lambda_k$ known (estimable from noise
  at the same station), the theorem holds with $Q_k$ replaced by $\Lambda_k^{1/2}Q_k$; for non-linear maps this is not
  proved here. With vertical transients and a random $Q$ the script recorded no violations (same table,
  row `heavy_pol`), but that is not a guarantee: a transient polarized along the subspace of a hypothesis
  would break it.
- **Independence of P and S** at the same station, as in AA.9.5.
- **A weaker alternative and why it is not enough.** If only per-station sign symmetry is assumed (the law of the
  record of station $i$ is invariant under $\mathbf y_i\mapsto-\mathbf y_i$, independent stations), the same argument
  with per-station Rademacher signs gives a valid bound without isotropy. But the conditional p-value of a sign
  test cannot go below $2^{1-n_s}$ (1/64 with 7 stations; the observed pattern is one of $2^{n_s-1}$
  equiprobable ones), which rules out the $p/N$ levels of a grid. With networks of 4 to 12 stations, which is the target
  regime, a sign randomization test is of no use; one by continuous rotations or by *time-slides* (§4.6 of the
  model design note) does not have that floor.

## AB.C.3 Absolute bounded influence implies saturation

**Proposition AB.C.3.** Let $\sigma:\mathbb R^{3n_b}\to[0,\infty)$ be any score (not necessarily radial) with finite
gross-error sensitivity,
$\mathrm{GES}=\sup_{\mathbf v,k,\mathbf w}|\sigma(\mathbf v)-\sigma(\mathbf v^{k\leftarrow\mathbf w})|<\infty$, where $\mathbf v^{k\leftarrow\mathbf w}$ replaces
block $k$ by $\mathbf w$. Then $\sigma(\mathbf v)\le\sigma(0)+n_b\,\mathrm{GES}$ for every $\mathbf v$: the score is bounded and, for
any event pattern $\mathbf w$, $\sigma(a\mathbf w)$ does not grow with $a$.

*Proof.* The blocks of $\mathbf v$ are replaced by zero one at a time; each step changes $\sigma$ by at most GES. $\square$

For the cap, $\mathrm{GES}=2c$ (AA.9.2) and the proposition recovers, with another constant, the saturation of AA.9.3. The
new content is that **no** score can have at the same time influence bounded uniformly in the rest of the
data and ordering by amplitude: the influence must be able to grow with the level of the other blocks. The two ways out
are a bound that grows slowly without being bounded (LOG6: influence $2c\sqrt{\log(1+M^2/c^2)}$, growth
$\sqrt{\log}$) or a bound **relative** to the background, which is what H5 of the model design note asks for and what
its contamination likelihood does (§4.4 of that note, whose properties (i) and (ii) are compatible with this
proposition precisely because its influence bound depends on $f^{\rm background}$).

## AB.C.4 ACAP: cap relative to the block median

**Definition.** Let $q(\mathbf v)$ be the $j$-th smallest of $\lVert\mathbf v_1\rVert,\dots,\lVert\mathbf v_{n_b}\rVert$, with $j=\lceil n_b/2\rceil$, and
$c_1>1$. ACAP projects each block onto the ball of radius $c_1q(\mathbf v)$:

$$\phi_A(\mathbf v)_k=\mathbf v_k\min\big(1,\ c_1q(\mathbf v)/\lVert\mathbf v_k\rVert\big),\qquad \sigma_A=\lVert Q_h^{\mathsf T}\phi_A(\mathbf v)\rVert .$$

**Proposition AB.C.4.**
(a) *Scale equivariance and non-saturation.* $\sigma_A(a\mathbf v)=|a|\sigma_A(\mathbf v)$, and $\sigma_A$ is $(1+c_1\sqrt{n_b})$-Lipschitz.
For an event pattern $\mathbf w$ with $\sigma_A(\mathbf w)>0$ over any background $\mathbf b$,
$\sigma_A(a\mathbf w+\mathbf b)\ge|a|\sigma_A(\mathbf w)-(1+c_1\sqrt{n_b})\lVert\mathbf b\rVert\to\infty$ linearly in $a$.
(b) *Bounded relative influence.* If $\mathbf v'$ differs from $\mathbf v$ in a set $K$ of $m\le\min(j-1,n_b-j)$ blocks, and $q_-$, $q_+$
are the $(j-m)$-th and the $j$-th smallest of the norms of the $n_b-m$ blocks outside $K$, then

$$|\sigma_A(\mathbf v)-\sigma_A(\mathbf v')|\le c_1\big[2\sqrt m\,q_++\sqrt{n_b-m}\,(q_+-q_-)\big],$$

finite whatever the substituted values. The breakdown point is $\min(j-1,n_b-j)/n_b$ (6 of 14 blocks).
(c) *Zero cost away from the extremes.* If $\max_k\lVert\mathbf v_k\rVert\le c_1q(\mathbf v)$, then $\phi_A(\mathbf v)=\mathbf v$ and $\sigma_A$ is the
uncapped score.
(d) *Tail.* ACAP is norm-adaptive, so Theorem AB.C.2 gives it its conditional threshold, with $\rho_k=\min(\lVert\mathbf v_k\rVert,c_1q)$.

*Proof.* (a) Norms, order statistic and projections scale with $|a|$, and the sign comes out of the norm. For
the Lipschitz constant: the projection onto a ball satisfies $\lVert P_{B(R)}\mathbf x-P_{B(R')}\mathbf x'\rVert\le\lVert\mathbf x-\mathbf x'\rVert+|R-R'|$, and
an order statistic of the norms is 1-Lipschitz with respect to $\max_k\lVert\mathbf v_k-\mathbf v'_k\rVert\le\lVert\mathbf v-\mathbf v'\rVert$; summing
over the blocks, $\lVert\phi_A(\mathbf v)-\phi_A(\mathbf v')\rVert\le\lVert\mathbf v-\mathbf v'\rVert+c_1\sqrt{n_b}\lVert\mathbf v-\mathbf v'\rVert$, and the contraction of Lemma
AA.9.1 carries it over to $\sigma_A$. The bound in $a$ is the Lipschitz inequality between $a\mathbf w+\mathbf b$ and $a\mathbf w$. (b) *Interlacing:* with
$c_{(1)}\le\dots\le c_{(n_b-m)}$ the clean norms, the $j$-th smallest of the full set lies in $[c_{(j-m)},c_{(j)}]$, because
among the $j$ smallest there are at most $m$ contaminated ones (hence at least $j-m$ clean ones below) and below $c_{(j)}$
there are at least $j$ values. Thus $q(\mathbf v),q(\mathbf v')\in[q_-,q_+]$. In the blocks of $K$, $\lVert\phi_A(\mathbf v)_k\rVert\le c_1q$ and likewise for $\mathbf v'$,
so that part contributes at most $\sqrt m\,c_1(q+q')\le2\sqrt mc_1q_+$; outside $K$ the blocks are equal and
the projection differs by at most $c_1|q-q'|\le c_1(q_+-q_-)$ per block. The contraction finishes the proof. (c) No block
exceeds the radius. (d) $\rho_k$ is a function of the norms. $\square$

**Reading.** ACAP behaves like the fixed cap when needed and like the identity when not:

- in Gaussian noise its effective radius is $c_1\bar q=7.2$ (with $c_1=3$), it almost never acts, and its null quantiles are
  those of the identity (table of AB.C.1);
- in the face of an isolated transient its influence is that of a cap at level $c_1q_+$, the local background, not that of its amplitude;
- in a saturated background (all norms at 20, above the fixed cap of 6) the cap stays pinned at its bound
  $c\sqrt{n_b}=22.45$ (20.6 with an event 100 times the background), LOG6 grows from 10 to 70 and ACAP grows linearly
  with the event, from 21 to 7,490, like the identity (`check_adaptive_out.txt`, C.4);
- the bound of (b) held in the 900 substitutions of the script (1, 3 and 6 blocks of amplitude up to $10^5$), with the
  maximum change at 26 % of the bound.

**What this does not resolve.** (1) The difference between ACAP and the fixed cap in R0 and C18 is a prediction, not a result:
the proposition says that ACAP does not saturate in a uniform background; if the background of a real sequence is very uneven across
stations, the median can fall on quiet stations and ACAP clips the noisy ones like a fixed cap. It is the
natural registrable test for phases 3 and 4. (2) The p-value of Corollary AB.C.2 is scale invariant and, for that reason,
**saturates** with strong events (it tends to a function of the coherence of the pattern). This is unavoidable by Proposition
AB.C.3 applied to a scale-invariant score, and it suggests separating roles: **deciding** with the conditional p-value
(calibrated against the background of that hypothesis) and **ranking** with $\sigma_A$ (which does not saturate). (3) Nothing in this part covers the
greedy deflation or the windows with events; that is Part B and the false-detection control of the complete detector (secondary
line 4 of the roadmap).

---

# Part D. Theory for the components of the new model

The model design note (§4) proposes replacing the cap
with a per-block contamination likelihood (§4.4), calibrating false detections with *time-slides* and adaptive conformal
inference (§4.6), and replacing the greedy deflation with a joint decomposition (§4.5). This part proves what can be
proved about each one and points out where the properties that the note marks as "inferred" need one more
assumption.

## AB.D.1 The per-block contamination likelihood

**Definition.** With the background density $f_k$ of block $k$ (estimated in the same noise regime), the
contamination rate $\varepsilon\in(0,1)$ and the width $s_k$ of the model error,

$$\ell_k(\mathbf m)=\log\frac{\varepsilon f_k(\mathbf v_k)+(1-\varepsilon)\,\mathcal N(\mathbf v_k;A_k\mathbf m,s_k^2I_3)}{f_k(\mathbf v_k)},\qquad
S_{\rm GLR}=\max_{\mathbf m}\sum_k\ell_k(\mathbf m),\qquad S_{\rm mix}=\log\int e^{\sum_k\ell_k(\mathbf m)}\,\pi(d\mathbf m),$$

with $\pi$ a prior distribution of the tensor fixed before seeing the data. It is the log-ratio of "the hypothesis
explains the block or the block is background" against "everything is background".

**Proposition AB.D.1.**
(a) *Influence.* For every $\mathbf m$, $\log\varepsilon\le\ell_k(\mathbf m)\le\log\big(\varepsilon+(1-\varepsilon)N^{\max}_k/f_k(\mathbf v_k)\big)$, with
$N^{\max}_k=(2\pi s_k^2)^{-3/2}$. If block $k$ is replaced by $\mathbf w$, $S_{\rm GLR}$ and $S_{\rm mix}$ change by at most
$\log\big(1+\tfrac{1-\varepsilon}{\varepsilon}N^{\max}_k/f_k(\mathbf w)\big)$ upwards (and the same with $\mathbf v_k$ downwards). The
influence is finite for each value but **not bounded**: it grows as $\log(1/f_k(\mathbf w))$, quadratically in
$\lVert\mathbf w\rVert$ if the background is modelled as Gaussian and as $(3+\nu)\log\lVert\mathbf w\rVert$ if it has polynomial tails of index
$\nu$, $f_k(\mathbf w)\sim C\lVert\mathbf w\rVert^{-(3+\nu)}$.
(b) *Non-saturation.* For an exactly predicted event, $\mathbf v_k=a\mathbf w_k$, with proportional model error
$s_k=s_0+\gamma a\lVert\mathbf w_k\rVert$ and a background with tails of index $\nu$, $\ell_k=\nu\log a+O(1)$ as $a\to\infty$; with fixed $s_k$,
$(3+\nu)\log a+O(1)$. If the background has scale $b_k$, $f_k(\mathbf v)=b_k^{-3}g(\mathbf v/b_k)$, the evidence depends on $a/b_k$: the
cut is **relative** to the background.

*Proof.* (a) The numerator lies between $\varepsilon f_k$ and $\varepsilon f_k+(1-\varepsilon)N^{\max}_k$. For $S_{\rm GLR}$, the change is at
most $\sup_{\mathbf m}[\ell_k(\mathbf w;\mathbf m)-\ell_k(\mathbf v_k;\mathbf m)]\le\sup\ell_k(\mathbf w;\cdot)-\log\varepsilon$; for $S_{\rm mix}$, the ratio of
the two integrals is bounded by the supremum of the ratio of the integrands, which is the same. The growth
rates are those of $\log(1/f_k)$. (b) $N^{\max}_k\sim(\gamma a\lVert\mathbf w_k\rVert)^{-3}(2\pi)^{-3/2}$ and
$1/f_k(a\mathbf w_k)\sim a^{3+\nu}\lVert\mathbf w_k\rVert^{3+\nu}/C$; the ratio grows as $a^\nu$. $\square$

**Reading.** (a) makes property (i) of §4.4 of the model design note precise: the influence bound is
$\log[(1-\varepsilon)/\varepsilon]$ **plus a term that grows with the amplitude of the transient** through the tail of the background. It is
compatible with Proposition AB.C.3 (nothing that does not saturate can have uniformly bounded influence). Concrete design
consequence: **the background must be modelled with polynomial tails**; a Gaussian background (or a mixture of Gaussians)
turns a transient that the hypothesis can fit into evidence quadratic in its amplitude, worse than the uncapped
score. (b) says that the tail index of the background relative to the model error decides non-saturation: with model
error proportional to the amplitude, the growth is $\nu\log a$, and if the background were as heavy as the model
error ($\nu\to0$) it would saturate again. With $\nu=3$ the script measures a slope of 3.00 over a dense background of scale
20 (`experiments/exploratory/THEORY-5/check_extensions_out.txt`, D.1), and the influence of a transient grows by 13.8 per decade of amplitude, the
$(3+\nu)\ln10$ of (a).

## AB.D.2 The mixture is an e-value

**Theorem AB.D.2.** If under $H_0$ the blocks are independent with densities exactly $f_k$, then for every
$\pi$ fixed in advance $\mathbb E_{H_0}[e^{S_{\rm mix}}]=1$, and by Markov's inequality $P_{H_0}(S_{\rm mix}\ge\log(1/p))\le p$.
The same holds if the integral is replaced by the average over $J$ draws $\mathbf m_j\sim\pi$ independent of the
data. If the true null density is $f_k^\ast$, $\mathbb E[e^{S_{\rm mix}}]=\prod_k\big(\varepsilon+(1-\varepsilon)\,\mathbb E_{\mathcal N(A_k\mathbf m,s_k^2)}[f_k^\ast/f_k]\big)$
averaged over $\pi$, so that a bound $f_k^\ast\le(1+\beta_k)f_k$ in the region where the model pulse has mass
gives an actual level of at most $p\prod_k(1+(1-\varepsilon)\beta_k)$.

*Proof.* For fixed $\mathbf m$, $\mathbb E_{f_k}[e^{\ell_k(\mathbf m)}]=\int(\varepsilon f_k+(1-\varepsilon)\mathcal N)\,d\mathbf v=1$; by independence the
product has mean 1; by Fubini (everything is non-negative) so do the mixture over $\pi$ and the Monte Carlo average. With
$f^\ast_k$ in place of $f_k$, the same computation gives the product written. $\square$

The threshold $\log(1/p)$ is **universal**: it does not depend on the network, the grid or the noise level, in exchange for the
background being well specified. Over $N$ hypotheses the union bound applies ($p/N$), and since averages of e-values are
e-values, the e-BH procedure (Wang and Ramdas, 2022) controls the proportion of false detections among the detections **without
assuming independence between hypotheses**, which is what §4.6 of the model design note asks of the conformal
p-values with BH (these need a positive-dependence condition that has not been checked on a grid). The
GLR ($\max$ over $\mathbf m$) is not an e-value and needs calibration.

*Verification* (`check_extensions_out.txt`, D.2; 14 blocks, $\varepsilon=0.05$, isotropic $t$ background with $\nu=3$ and lognormal
scales, $\pi=\mathcal N(0,10^2I)$, 4000 draws, 4000 null replicates). With the background well specified there was
no exceedance at $p=0.1$, 0.01 and 0.001, nor with an actual background twice as wide. With an actual background with heavier
tails ($\nu=1.5$) the mixture still had no exceedances and the GLR exceeded at 0.0013, 0.0008 and 0.0003 (below $p$ here,
but not guaranteed). The other side: with a diffuse prior the e-value is very conservative (the empirical mean of
$e^{S_{\rm mix}}$ is almost 0 because its mean of 1 is supplied by rare, enormous values); its power depends on $\pi$, and a $\pi$
centred on the tensors of the area (for example, that of the calibration events) is the natural choice.

## AB.D.3 *Time-slides*: which null they calibrate

**Proposition AB.D.3.** Let the record of a circular window of length $T$ samples be $\mathbf y=(\mathbf y_1,\dots,\mathbf y_{n_s})$, and
$G=\mathbb Z_T^{n_s}$ the group of independent circular shifts per station. If under $H_0$ the records of the
stations are independent and each one is stationary in the window, then $g\mathbf y\overset{d}{=}\mathbf y$ for every $g\in G$, and
for any statistic $M$ (for example, the maximum over the grid of the complete detector, deflation included) the
p-value $\hat p=(1+\#\{b:M(g_b\mathbf y)\ge M(\mathbf y)\})/(B+1)$, with $g_b$ uniform and independent in $G$, satisfies
$P_{H_0}(\hat p\le u)\le u$.

*Proof.* With $g_0$ uniform on $G$ independent of everything, $\mathbf y\overset{d}{=}g_0\mathbf y$, and $(g_0\mathbf y,g_1\mathbf y,\dots,g_B\mathbf y)=(g_0\mathbf y,(g_1g_0^{-1})g_0\mathbf y,\dots)$
where the $g_bg_0^{-1}$ are uniform and independent of each other and of $g_0\mathbf y$; the sequence $(M(\mathbf y),M(g_1\mathbf y),\dots)$ is
exchangeable, and the rank of $M(\mathbf y)$ is uniform up to ties, which only make the p-value conservative
(Hemerik and Goeman, 2018). $\square$

**What they calibrate and what they do not.** The null is **global**: no coherent source in the window. With a real event in
the window the stations are no longer independent and the null is false; therefore the *time-slides* do not control the
false detections that **a real event causes** in wrong hypotheses (lobes, duplicates, the residual of
Proposition AB.B.6), which are coherent by construction and which a shift destroys. They do control the false detections
caused by noise, transients and the incoherent energy of the background, in the noise regime of that same hour. In an
aftershock sequence this means: the *time-slides* calibrate the first step of the detector (is there anything coherent?),
and the false detections inside windows with events need another tool (Part B, the FROC, or a conditional null
that deflates the events already accepted before shifting). Unlike the sign test of §AB.C.2 there is no practical
floor on the p-value ($|G|=T^{n_s}$), and shifts larger than the maximum *moveout* are not needed for validity
under $H_0$; with events present, small shifts retain part of their coherence and make the p-value
conservative, not liberal. For adaptive conformal inference (Gibbs and Candès, 2021) the guarantee is deterministic and long-run:
with $\alpha_{t+1}=\alpha_t+\gamma(\alpha-\mathrm{err}_t)$, $|\tfrac1T\sum_t\mathrm{err}_t-\alpha|\le(\max(\alpha_1,1-\alpha_1)+\gamma)/(\gamma T)$, without
exchangeability, but marginal in time and not per window.

## AB.D.4 Joint decomposition: what certificates can and cannot promise

**Proposition AB.D.4.** (a) For two hypotheses with orthonormal bases $Q_1,Q_2$ and $\Gamma=[Q_1\ Q_2]$,
$\sigma_{\min}(\Gamma)^2=1-\cos\theta_1$, with $\theta_1$ the smallest principal angle between their ranges. (b) For two
co-located hypotheses with origins separated by $\Delta\tau$ and orthogonal P and S groups, $1-\cos\theta_1\le1-a(\Delta\tau)=\omega^2\Delta\tau^2/2+O(\Delta\tau^4)$, so $\sigma_{\min}(\Gamma)\le\omega\Delta\tau/\sqrt2\,(1+o(1))$, and the
sufficient conditions of AA.6(a) ($\lVert\mathbf w\rVert\le\zeta\lambda/2$, $\lambda<\zeta\sigma^3a_{\min}/(6\sqrt s)$) require
$a_{\min}>12\sqrt2\,\lVert\mathbf w\rVert/(\zeta^2\sigma^3)$, that is, an SNR of order $\Delta\tau^{-3}\zeta^{-2}$, against the necessary $\Delta\tau^{-1}$ of
Corollary AB.A.7.

*Proof.* (a) $\Gamma^{\mathsf T}\Gamma=\big(\begin{smallmatrix}I&C\\C^{\mathsf T}&I\end{smallmatrix}\big)$ with $C=Q_1^{\mathsf T}Q_2$, whose singular values are
the principal cosines $\cos\theta_i$; its eigenvalues are $1\pm\cos\theta_i$. (b) $\cos\theta_1=\max\langle\mathbf u_1,\mathbf u_2\rangle$ over unit vectors of the two ranges. Taking the same
tensor in both hypotheses, each arrival is the same pulse shifted by $\Delta\tau$, and with orthogonal groups
$\langle\mathbf u_1,\mathbf u_2\rangle=a(\Delta\tau)$; the expansion of $a$ gives the bound. Substituting
$s=2$ and $\lambda\ge2\lVert\mathbf w\rVert/\zeta$ into AA.6(a) yields the condition. $\square$

The script measures $1-\cos\theta_1=8.88\times10^{-4},\ 9.85\times10^{-3},\ 8.70\times10^{-2}$ for $\Delta\tau=3,10,30$ ms, slope 2.00
(`check_extensions_out.txt`, D.4); at 3 ms it coincides with $\omega^2\Delta\tau^2/2=8.88\times10^{-4}$, so the bound of (b) is
tight: the best aligned pair is the same tensor shifted. **Consequence:** no certificate of coherence or pre-certificate type, for the
group BLASSO, BOMP or a joint decomposition with the same dictionary, can guarantee the separation of two
overlapping events near the information limit; the gap is a factor $\Delta\tau^{-2}$ in the SNR. This is the
known phenomenon of super-resolution (convex guarantees require minimum separation), here with its exponent. The
ADCG-type joint decomposition (§4.5 of the model design note) is the right algorithm for approaching the
limit, because it re-optimizes continuous positions and times, but its guarantees in this regime remain in the state of
Conjecture AA.1. Two things can be stated: (i) a subtraction with refitted sub-sample time and amplitude (the minimum
of a joint re-optimization) removes the term of Proposition AB.B.6, so that the residual of a found event
is only the projection of the noise and of the model error; (ii) with a dictionary of empirical Green's functions, the
mismatch $\eta$ of an atom against the signal enters Lemma AB.B.2 as $-\eta\lVert\mathbf Y\rVert$ and does not enter Theorem AB.A.3
(the timing-tolerant class does not depend on the dictionary), which cleanly separates what the dictionary changes
(detection) from what it does not (identifiability).

**Conjecture AB.D.5 (target statement of §6 of the model design note, reformulated).** For two events with
$K_{\rm mec}>0$, the joint decomposition with joint reading selection and contamination likelihood separates the
two events with probability $\ge1-e$ when $\delta_{\rm TT}\ge C\,z_e$ for a constant $C$ that does not depend on $\Delta\tau$. It is the
claim that the algorithm attains the information limit up to constants; Proposition AB.D.4 says that it cannot
be proved with coherence certificates, and the candidate tool is a local analysis of the joint
re-optimization in the neighbourhood of the truth (like that of Duval and Peyré, 2015, for the scalar BLASSO), with $F_{\rm pol}$ of
Corollary AB.A.6 as the curvature.

---

## AB.5 Integration into the repository

The integration was carried out as follows:

1. **Thesis Section 4 and AA.9.5(b):** both now cite Proposition AB.C.1, with $\kappa$ replaced by 1 in the bound; at
   $p=10^{-4}$, 28.50 became 14.59, and "conservative by a factor of four" became "by a factor of about two".
   `experiments/exploratory/THEORY-5/tail_check_c1.py` prints both bounds (the frozen
   `code/theory_tail_check.py` and its output are unchanged).
2. **Appendix AB:** this appendix. Part A is the primary line of phase 5, and its gate (complete statements,
   proofs and numerical verification) is met for AB.A.1 to AB.A.7; AB.A.8 remains a conjecture.
3. **Error register:** Proposition AB.B.6 (residual of the deflation by subtraction at the sample read) was not among
   E1 to E56; E52 deals with the rounding of the exact-time readers, not with the deflation. Proposition AB.B.6 is
   proposed as entry E62 of the error register, with the figure measured in §AB.B.5; the entry and a sub-sample
   deflation convention are under review in pull request #5.
4. **Atlas:** splitting the < 0.5 s cell and recording $K_{\rm mec}/E$ per pair (§AB.A.7) is pending in the atlas card
   (Phase 3).
5. **New model (model design note):** the following consequences are kept as the design inputs for Phase 4: model the
   background of the contamination likelihood with polynomial tails (Proposition AB.D.1); use the mixture version as an
   e-value with e-BH for the detections (Theorem AB.D.2); use the *time-slides* only for the global null and treat the
   false detections in windows with events separately (Proposition AB.D.3); in the joint decomposition, subtract with
   sub-sample time (Proposition AB.B.6) and do not expect coherence certificates in the overlapping regime
   (Proposition AB.D.4).

## References

- Archimedes' hat-box theorem: the projection of the uniform measure on $S^2$ onto an axis is uniform.
- Boucheron, S., Lugosi, G., Massart, P. (2013). *Concentration Inequalities*. Oxford University Press (Gaussian
  concentration of Lipschitz functions, maximum of sub-Gaussians).
- Hsu, D., Kakade, S. M., Zhang, T. (2012). A tail inequality for quadratic forms of subgaussian random vectors.
  *Electronic Communications in Probability* 17, 52.
- Le Cam, L. (1973). Convergence of estimates under dimensionality restrictions. *Annals of Statistics* 1, 38–53
  (two-point bound via total variation).
- Hampel, F. R. (1974). The influence curve and its role in robust estimation. *JASA* 69, 383–393 (gross-error
  sensitivity, breakdown point).
- Eldar, Y. C., Kuppinger, P., Bölcskei, H. (2010). Block-sparse signals: uncertainty relations and efficient recovery.
  *IEEE Trans. Signal Processing* 58, 3042–3054 (block coherence for BOMP; Lemma AB.B.3 is its timing-tolerant
  analogue with peak reading).
- Candès, E., Fernandez-Granda, C. (2014). Towards a mathematical theory of super-resolution. *CPAM* 67, 906–956 (the
  scaling laws of the minimum resolvable separation with the SNR).
- Wang, R., Ramdas, A. (2022). False discovery rate control with e-values. *JRSS B* 84, 822–852.
- Hemerik, J., Goeman, J. (2018). Exact testing with random permutations. *TEST* 27, 811–825.
- Gibbs, I., Candès, E. (2021). Adaptive conformal inference under distribution shift. *NeurIPS*.
- Duval, V., Peyré, G. (2015). Exact support recovery for sparse spikes deconvolution. *Found. Comput. Math.* 15, 1315–1355.
- Appendix AA of the supplementary evidence and the model design note (state of the art and design of the new model)
  of this project.
