# C8BC-Supplementary
C8BC Supplementary Material
Cyclic Eight-Brane Collision Cosmology: A CPT-Symmetric Multi-Brane Framework Unifying the CMB Cold Spot, Baryon Asymmetry, Dark Matter, and Dark Energy
Author: Aaron Riley Hurst
Status: Submitted to the Journal of Cosmology and Astroparticle Physics (JCAP), 2026
Companion paper: Hurst, A.R. (2026). The Rapid Planetary Disassembly Hypothesis: 22 Lines of Evidence for a Catastrophic Disruption at 2.19 AU 4.568 Billion Years Ago. Under editorial review, Icarus. Monte Carlo code: https://github.com/Aaronius0000/Simulated-Asteroid-Trajectories-from-RPD-Explosion-Events

Overview
This repository contains supplementary code, derivation notes, and search templates accompanying the C8BC paper. The Cyclic Eight-Brane Collision Cosmology (C8BC) is a theoretical framework in which the Big Bang originated in the collision of two CPT-mirror gravitational singularities from a previous cosmic cycle, potentially producing eight stable brane configurations organized by the rotational symmetry of the compact extra-dimensional manifold. The framework simultaneously addresses five major unsolved problems in cosmology — the CMB Cold Spot and Eridanus Supervoid, the baryon asymmetry, dark matter distribution and granularity, the dark energy fine-tuning problem and equation of state evolution, and the Hubble tension — through a single physical mechanism, and makes sixteen specific falsifiable predictions testable with instruments currently operational or under construction.
The framework additionally predicts the existence of other CPT-mirror cosmological structures in the shared extra-dimensional bulk, detectable as anomalous far-field dark gravity contributions from directions unrelated to the Cold Spot collision axis — Prediction 16, the Grand Cosmic Horizon.

Repository Contents
Code S2 — Northern Cold Spot Detection Probability
File: code_S2_northern_cold_spot_detection.py
A Python Monte Carlo simulation estimating the probability that CMB-S4 detects the C8BC predicted Northern Cold Spot as a function of signal amplitude and instrument noise, using a matched spherical cap filter on simulated CMB sky realizations.
The C8BC framework predicts two distinct CMB temperature deficit features at antipodal sky positions — the Southern Cold Spot (the observed Eridanus Supervoid bilateral collision boundary scar) and a structurally distinct Northern Cold Spot (a leading-edge expansion stress feature at the antipodal position). The Northern Cold Spot is predicted with amplitude 50 to 140 microkelvin, a shallower and more diffuse internal profile than the Southern Cold Spot, and directional E-mode polarization asymmetry along the expansion axis rather than the bilateral symmetry of the Southern Cold Spot.
The code computes:

Detection probability as a function of signal amplitude across the C8BC predicted range (50 to 140 microkelvin) at CMB-S4 noise levels
Null hypothesis false positive rate at the 3-sigma detection threshold
Structural asymmetry discrimination probability — the probability of correctly identifying the C8BC asymmetric model over a symmetric bubble collision alternative given both features are detected

Dependencies: numpy, scipy, matplotlib (all standard, no external data required)
Usage:
bashpython code_S2_northern_cold_spot_detection.py
Outputs detection probability curve, signal vs null comparison, and structural asymmetry discrimination figure saved as c8bc_northern_cold_spot_detection.png.

Template S6 — SKA RPD Gravitational Wave Search Template
File: template_S6_SKA_RPD_search.pdf
A preliminary search template for Square Kilometre Array pulsar timing array analysis teams targeting the coordinated multi-brane Rapid Planetary Disassembly gravitational wave signature of Predictions 7A and 7B.
The C8BC framework predicts two distinct gravitational wave signals from the coordinated multi-brane RPD Terra-Forging event 4.568 billion years ago — Signal 7A from the Wing One RPD at cosmological lookback time 4.568 billion years arriving from the galactic plane direction, and Signal 7B from the Wing One Orthogonal bulk-propagating RPD component arriving from the Cold Spot direction at cosmological lookback time 4.568 billion minus 2.77 million years. The phase offset between the two signals directly measures the extra-dimensional bulk separation distance d_bulk, providing one of three independent measurements of this parameter predicted by the C8BC framework.
Template specifies: Expected signal frequency (~6.7 nanohertz), amplitude order-of-magnitude estimates, sky direction coordinates for Signal 7A (galactic plane) and Signal 7B (Cold Spot direction, approximately RA 03h 15m, Dec −19°), predicted lookback time separation (~2.77 million years), and a five-step recommended SKA analysis strategy for near-field/far-field d_bulk directional decomposition.
Note: This template is explicitly preliminary. Signal amplitude estimates will be refined when the formal parameter fitting pipeline (Code S3, in preparation) is complete and when the RPD companion paper's Monte Carlo constraints on the disruption event energy are finalized. The template should not be used as a formal SKA search specification without the refined amplitude calculations.

Supplementary Note S5 — Dark Solar System Gravitational Anomaly Calculations
File: note_S5_dark_solar_system_calculations.pdf
Detailed numerical calculations for Prediction 10 Components A, B, and C — the three-component dark gravity structure predicted in the solar system from natural bodies, engineered bodies, and the Dark Belt distributed profile.
Contains:

Component A: Dark Sun gravitational contribution at the Earth-Sun distance (ρ_DM ≈ 2.46 × 10⁻¹⁰ kg m⁻³ at 1 AU, 14 orders of magnitude above smooth halo background), Dark Jupiter, Dark Saturn contributions
Component B Scenario A: Dark Earth GM excess calculation (ΔGM/GM ≈ 5.2 × 10⁻⁶ above Component A background), Dark Nibiru mass discrepancy (ΔM/M ≈ 5 × 10⁻⁶ of Planet Nine's mass)
Component B Scenario B: null prediction and discriminant conditions
Component C: Dark Belt distributed dark gravity profile spanning 2.0–3.5 AU with peak near 2.19 AU and 0.1–0.3 AU positional offsets from current large belt member positions

Important correction noted: The Planet Nine mass discrepancy prediction in the main text of the companion paper stated 10 to 20 percent. The derivation-consistent value computed in this note is approximately 5 × 10⁻⁶ of Planet Nine's mass — four orders of magnitude smaller. The main text has been corrected accordingly.

Supplementary Note S1 — Monte Carlo Posterior Analysis Description
File: note_S1_monte_carlo_description.pdf
Description of the Monte Carlo parameter constraint methodology used to establish the order-of-magnitude parameter estimates in Table 1 of the main paper. Parameter values in Table 1 are order-of-magnitude estimates constrained by the combined Planck 2018, DESI DR2, and SH0ES observational bounds and are not formal best-fit results from a likelihood analysis. The formal parameter fitting pipeline is identified as future work (Code S3, in preparation).

In Preparation
The following supplementary items are identified as priority future work and will be deposited to this repository upon completion:
Code S3 — Parameter Fitting Pipeline
Joint likelihood analysis of C8BC framework predictions against the full Planck 2018 CMB power spectrum, DESI DR2 BAO distance measurements, and SH0ES H₀ data over the six-dimensional C8BC parameter space. Will produce formal best-fit values and posterior uncertainties for all six parameters replacing the current order-of-magnitude estimates in Table 1.
Figures S4 — Corner Plots and Posterior Distributions
Posterior probability distributions and two-dimensional joint posteriors for all six C8BC parameters from Code S3. Downstream of Code S3 and will be deposited simultaneously.
Readers with expertise in Markov Chain Monte Carlo cosmological parameter estimation who are interested in contributing to the Code S3 pipeline are invited to contact the author.

The Five Anomalies and Sixteen Predictions
The C8BC framework addresses five major unsolved problems in cosmology simultaneously:
AnomalyC8BC MechanismCMB Cold Spot and Eridanus SupervoidBilateral collision boundary scar from matter evacuated into both expanding wings simultaneouslyBaryon asymmetryNatural consequence of CPT collision geometry — each wing inherits opposite charge asymmetryDark matter distribution and granularityThree-source dark gravity from neighboring matter-wing branes with proximity-ordered coupling strengths 1.00 : 0.60 : 1.00Dark energy fine-tuning and equation of stateEight-brane vacuum cancellation residual plus bulk tension from expanding wing separationHubble tensionComponent Two bulk tension modifying integrated expansion history; Component Three KK correction contributing directional H₀ anisotropy
The framework makes sixteen specific falsifiable predictions. All fifteen directional predictions share alignment with the Cold Spot collision axis. Prediction 16 tests for far-field dark gravity contributions from uncorrelated directions as evidence for other CPT-mirror cosmological structures in the shared extra-dimensional bulk.
Primary falsifiers — any one rules out the C8BC:

CMB-S4 finds no Northern Cold Spot at the antipodal sky position at greater than 3-sigma (early 2030s)
Euclid finds no preferred axis in large-scale structure aligned with the Cold Spot at greater than 3-sigma
DESI confirms w = −1 exactly with no time variation at greater than 3-sigma
Euclid finds no sub-kiloparsec dark gravity-baryon cross-correlation at greater than 3-sigma
Multiple directional signals detected with mutually inconsistent preferred axes


Brane Architecture Summary
The C8BC framework proposes up to eight stable brane configurations at equal 45-degree intervals in the compact extra-dimensional spin space, plus a spin-neutral center position at the geometric midpoint equidistant from all eight:
PositionNameRole0°Wing OneOur universe — reference frame+45°Plus-45-VerseDominant anisotropic dark gravity source (coupling 1.00)+90°Wing One OrthogonalMinor smooth isotropic dark gravity source (coupling 0.60)−45° / 315°Minus-45-VerseDominant anisotropic dark gravity source (coupling 1.00)180°Wing TwoCPT-mirror antimatter twin+225°Plus-45-Anti-VerseWing Two's closest neighbor+270°Wing Two OrthogonalWing Two's smooth dark gravity source+135°Minus-45-Anti-VerseMatter-antimatter boundary bridge braneCenterSpin-Neutral CenterEquidistant from all eight branes — the ninth position
The matter wing spans strictly less than 180 degrees of spin space — approached asymptotically as the cyclic brane accumulation adds finer subdivisions toward the antimatter boundary with each successful generational reset, each requiring external matter acquisition from neighboring dying universes in the surrounding bulk.

Cyclic Generational History
GenerationConfigurationDark MatterHalo AlignmentΛ ResidualFirstTwo-wing, two-braneNoneNone~UnmitigatedSecondTwo-wing, four-braneSingle isotropic sourceNoneReducedThird (this universe)Two-wing, eight-braneThree-source, anisotropic dominantPresentDramatically reducedBeyondTwo-wing, beyond eight-braneIncreasingly fine multi-sourceStrengthenedAsymptotically approaching zero
Each generational advancement requires external matter from neighboring dying universes acquired by the spin-neutral center occupants between resets. The universe does not automatically improve — it improves only if sufficient external matter is acquired.

Companion Repository
The Monte Carlo statistical analysis for the Rapid Planetary Disassembly companion paper (under review at Icarus) is maintained separately at:
https://github.com/Aaronius0000/Simulated-Asteroid-Trajectories-from-RPD-Explosion-Events
The RPD and C8BC frameworks are independently testable but mutually reinforcing. The RPD event parameters constrain the expected amplitude and frequency of the coordinated gravitational wave signal in C8BC Predictions 7A and 7B. A joint calculation connecting the two frameworks for the SKA search template is identified as a priority future work item.

Citation
If you use this code or supplementary material, please cite:
Hurst, A.R. (2026). Cyclic Eight-Brane Collision Cosmology: A CPT-Symmetric Multi-Brane Framework Unifying the CMB Cold Spot, Baryon Asymmetry, Dark Matter, and Dark Energy. Submitted to JCAP.
And this repository:
Hurst, A.R. (2026). C8BC Supplementary Material. GitHub. https://github.com/Aaronius0000/C8BC-Supplementary

Contact
Aaron Riley Hurst
Independent researcher, Brandon, Florida
Correspondence regarding the C8BC framework, collaboration inquiries regarding the Code S3 parameter fitting pipeline, and requests for the full manuscript preprint are welcome.

License
Code released under MIT License. Supplementary notes released under CC BY 4.0. Both licenses permit free use with attribution.

"We will know within a decade whether the universe left a forwarding address."
