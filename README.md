"""
Supplementary Code S2
Northern Cold Spot Detection Probability Under C8BC versus Null Hypothesis

Estimates the probability that CMB-S4 detects the Northern Cold Spot
predicted by the C8BC framework as a function of signal amplitude
and instrument noise, using a matched spherical cap filter on a
simulated CMB sky.

Dependencies: numpy, scipy, matplotlib
No external data required — all inputs are parameterized.

Author: A. R. Hurst
Paper: Cyclic Eight-Brane Collision Cosmology (C8BC)
Date: 2026
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Instrument and sky parameters ────────────────────────────────────────────

# CMB-S4 noise level in microkelvin per square-degree pixel
# Reference: CMB-S4 Science Book arXiv:1610.02743
NOISE_UK_PER_SQRDEG = 1.0          # CMB-S4 target noise level

# Sky pixel resolution in degrees
PIXEL_SIZE_DEG = 0.5

# Angular radius of the Cold Spot feature in degrees
# Southern Cold Spot observed radius approximately 5 degrees
COLD_SPOT_RADIUS_DEG = 5.0

# Number of sky pixels in the analysis region
# We work on a local patch around the predicted Northern Cold Spot position
# rather than a full-sky simulation for computational tractability
PATCH_SIZE_DEG = 30.0              # Analysis patch radius in degrees
N_PIXELS_1D = int(2 * PATCH_SIZE_DEG / PIXEL_SIZE_DEG)
N_PIXELS_TOTAL = N_PIXELS_1D ** 2

# ── C8BC signal parameters ────────────────────────────────────────────────────

# C8BC predicted Northern Cold Spot amplitude range in microkelvin
# Southern Cold Spot observed: 70-140 microkelvin deficit
# Northern Cold Spot predicted: 50-140 microkelvin (shallower, more diffuse)
SIGNAL_AMPLITUDE_MIN_UK = 50.0
SIGNAL_AMPLITUDE_MAX_UK = 140.0
SIGNAL_AMPLITUDE_CENTRAL_UK = 80.0    # Central estimate for detection calc

# Northern Cold Spot is predicted to be more diffuse than Southern
# Profile shallower by approximately 30-50 percent
NORTHERN_DIFFUSE_FACTOR = 0.7         # 70 percent of Southern gradient steepness

# ── Simulation parameters ─────────────────────────────────────────────────────

N_MONTE_CARLO = 10000              # Number of Monte Carlo realizations
DETECTION_THRESHOLD_SIGMA = 3.0    # Detection threshold in sigma
RANDOM_SEED = 42

np.random.seed(RANDOM_SEED)


# ── Helper functions ──────────────────────────────────────────────────────────

def gaussian_cap_profile(pixel_distances_deg, cap_radius_deg,
                          amplitude_uk, diffuse_factor=1.0):
    """
    Generate a spherical cap temperature deficit profile.

    Parameters
    ----------
    pixel_distances_deg : array
        Angular distance of each pixel from cap center in degrees
    cap_radius_deg : float
        Angular radius of the cap feature in degrees
    amplitude_uk : float
        Central temperature deficit in microkelvin (positive = deficit)
    diffuse_factor : float
        Controls boundary gradient steepness.
        1.0 = Southern Cold Spot steepness
        < 1.0 = Northern Cold Spot more diffuse boundary

    Returns
    -------
    array : temperature deficit in microkelvin at each pixel
    """
    # Gaussian profile with width set by cap radius and diffuse factor
    # Steeper boundary for Southern (diffuse_factor=1.0)
    # Shallower boundary for Northern (diffuse_factor < 1.0)
    sigma = cap_radius_deg / (2.0 * diffuse_factor)
    profile = amplitude_uk * np.exp(
        -0.5 * (pixel_distances_deg / sigma) ** 2
    )
    return profile


def simulate_cmb_patch(n_pixels_1d, pixel_size_deg, noise_uk_per_sqrdeg):
    """
    Simulate a CMB temperature patch with instrument noise.
    Uses a simple white noise model — adequate for detection threshold
    estimation. A full simulation would include the CMB primary anisotropy
    power spectrum but this does not significantly affect the
    matched-filter detection threshold at the Cold Spot angular scale.

    Parameters
    ----------
    n_pixels_1d : int
        Number of pixels along one side of the square patch
    pixel_size_deg : float
        Pixel size in degrees
    noise_uk_per_sqrdeg : float
        Instrument noise in microkelvin per square-degree pixel

    Returns
    -------
    array : 2D array of simulated CMB temperatures in microkelvin
    """
    noise_per_pixel = noise_uk_per_sqrdeg * pixel_size_deg
    return np.random.normal(0, noise_per_pixel, (n_pixels_1d, n_pixels_1d))


def build_pixel_distance_map(n_pixels_1d, pixel_size_deg):
    """
    Build a map of angular distances from patch center for each pixel.

    Returns
    -------
    array : 2D array of angular distances in degrees from patch center
    """
    center = n_pixels_1d / 2.0
    x_indices = np.arange(n_pixels_1d)
    y_indices = np.arange(n_pixels_1d)
    xx, yy = np.meshgrid(x_indices, y_indices)
    distances = np.sqrt((xx - center)**2 + (yy - center)**2) * pixel_size_deg
    return distances


def matched_filter_statistic(sky_map, signal_template):
    """
    Compute the matched filter detection statistic.

    The matched filter maximizes signal-to-noise for a known signal
    template in Gaussian noise. The statistic is the normalized
    cross-correlation of the sky map with the signal template.

    Parameters
    ----------
    sky_map : 2D array
        Observed or simulated sky temperature map
    signal_template : 2D array
        Expected signal profile (normalized)

    Returns
    -------
    float : detection statistic in units of noise sigma
    """
    template_normalized = signal_template / np.sqrt(
        np.sum(signal_template ** 2)
    )
    statistic = np.sum(sky_map * template_normalized)
    # Normalize by noise level per pixel
    noise_per_pixel = np.std(sky_map[signal_template < 0.01 *
                                      np.max(signal_template)])
    if noise_per_pixel == 0:
        noise_per_pixel = np.std(sky_map) + 1e-10
    return statistic / (noise_per_pixel *
                        np.sqrt(np.sum(template_normalized ** 2)))


# ── Main detection probability calculation ────────────────────────────────────

def compute_detection_probability(signal_amplitude_uk,
                                   diffuse_factor=NORTHERN_DIFFUSE_FACTOR,
                                   n_realizations=N_MONTE_CARLO,
                                   threshold_sigma=DETECTION_THRESHOLD_SIGMA,
                                   verbose=False):
    """
    Estimate the probability that CMB-S4 detects the Northern Cold Spot
    at a given signal amplitude.

    Parameters
    ----------
    signal_amplitude_uk : float
        True Northern Cold Spot amplitude in microkelvin
    diffuse_factor : float
        Boundary gradient factor relative to Southern Cold Spot
    n_realizations : int
        Number of Monte Carlo sky realizations
    threshold_sigma : float
        Detection threshold in sigma
    verbose : bool
        Print progress

    Returns
    -------
    float : detection probability (fraction of realizations above threshold)
    float : uncertainty on detection probability (binomial standard error)
    """
    distance_map = build_pixel_distance_map(N_PIXELS_1D, PIXEL_SIZE_DEG)

    # Build the signal template — the expected Northern Cold Spot profile
    signal_template = gaussian_cap_profile(
        distance_map,
        COLD_SPOT_RADIUS_DEG,
        signal_amplitude_uk,
        diffuse_factor
    )

    detections = 0

    for i in range(n_realizations):
        if verbose and i % 1000 == 0:
            print(f"  Realization {i}/{n_realizations}")

        # Simulate noise-only background
        noise_map = simulate_cmb_patch(
            N_PIXELS_1D, PIXEL_SIZE_DEG, NOISE_UK_PER_SQRDEG
        )

        # Add signal to noise
        sky_with_signal = noise_map - signal_template  # deficit = negative

        # Apply matched filter
        detection_stat = matched_filter_statistic(
            sky_with_signal, signal_template
        )

        if detection_stat >= threshold_sigma:
            detections += 1

    prob = detections / n_realizations
    # Binomial standard error
    uncertainty = np.sqrt(prob * (1 - prob) / n_realizations)

    return prob, uncertainty


def compute_null_false_positive_rate(n_realizations=N_MONTE_CARLO,
                                      threshold_sigma=DETECTION_THRESHOLD_SIGMA):
    """
    Estimate the false positive rate under the null hypothesis —
    the probability of claiming a detection when no signal is present.

    Returns
    -------
    float : false positive rate
    """
    distance_map = build_pixel_distance_map(N_PIXELS_1D, PIXEL_SIZE_DEG)

    # Template is the Northern Cold Spot profile shape with unit amplitude
    # The filter searches for this shape regardless of amplitude
    signal_template = gaussian_cap_profile(
        distance_map,
        COLD_SPOT_RADIUS_DEG,
        amplitude_uk=1.0,
        diffuse_factor=NORTHERN_DIFFUSE_FACTOR
    )

    false_positives = 0

    for _ in range(n_realizations):
        # Noise only — no signal
        noise_map = simulate_cmb_patch(
            N_PIXELS_1D, PIXEL_SIZE_DEG, NOISE_UK_PER_SQRDEG
        )
        detection_stat = matched_filter_statistic(noise_map, signal_template)
        if detection_stat >= threshold_sigma:
            false_positives += 1

    return false_positives / n_realizations


# ── Scan over signal amplitudes ───────────────────────────────────────────────

def amplitude_scan(amplitude_range_uk=None, n_points=20):
    """
    Compute detection probability as a function of signal amplitude
    across the C8BC predicted range.

    Parameters
    ----------
    amplitude_range_uk : tuple or None
        (min, max) amplitude range in microkelvin.
        Defaults to C8BC predicted range (50, 140).
    n_points : int
        Number of amplitude values to evaluate

    Returns
    -------
    amplitudes : array of amplitude values
    probabilities : array of detection probabilities
    uncertainties : array of probability uncertainties
    """
    if amplitude_range_uk is None:
        amplitude_range_uk = (SIGNAL_AMPLITUDE_MIN_UK,
                               SIGNAL_AMPLITUDE_MAX_UK)

    amplitudes = np.linspace(amplitude_range_uk[0],
                              amplitude_range_uk[1], n_points)
    probabilities = np.zeros(n_points)
    uncertainties = np.zeros(n_points)

    print(f"Running amplitude scan over {n_points} values "
          f"from {amplitude_range_uk[0]} to {amplitude_range_uk[1]} uK...")

    for i, amp in enumerate(amplitudes):
        prob, unc = compute_detection_probability(amp)
        probabilities[i] = prob
        uncertainties[i] = unc
        print(f"  Amplitude {amp:.1f} uK: "
              f"P(detect) = {prob:.3f} +/- {unc:.3f}")

    return amplitudes, probabilities, uncertainties


# ── Structural asymmetry test ─────────────────────────────────────────────────

def structural_asymmetry_discriminant(n_realizations=N_MONTE_CARLO):
    """
    Estimate the probability of distinguishing the C8BC structural asymmetry
    prediction from a symmetric bubble collision model.

    The C8BC predicts the Northern Cold Spot is more diffuse than the
    Southern (diffuse_factor ~ 0.7). A symmetric model predicts identical
    profiles (diffuse_factor = 1.0).

    Returns
    -------
    float : probability of correctly identifying C8BC asymmetry over
            symmetric alternative given both features are detected
    """
    distance_map = build_pixel_distance_map(N_PIXELS_1D, PIXEL_SIZE_DEG)

    # True C8BC signal — more diffuse Northern Cold Spot
    true_signal = gaussian_cap_profile(
        distance_map,
        COLD_SPOT_RADIUS_DEG,
        SIGNAL_AMPLITUDE_CENTRAL_UK,
        diffuse_factor=NORTHERN_DIFFUSE_FACTOR    # 0.7 — more diffuse
    )

    # Symmetric alternative template — same profile as Southern Cold Spot
    symmetric_template = gaussian_cap_profile(
        distance_map,
        COLD_SPOT_RADIUS_DEG,
        SIGNAL_AMPLITUDE_CENTRAL_UK,
        diffuse_factor=1.0                         # symmetric
    )

    correct_identifications = 0

    for _ in range(n_realizations):
        noise_map = simulate_cmb_patch(
            N_PIXELS_1D, PIXEL_SIZE_DEG, NOISE_UK_PER_SQRDEG
        )
        observed = noise_map - true_signal

        # Compute match to each template
        stat_c8bc = matched_filter_statistic(observed, true_signal)
        stat_symmetric = matched_filter_statistic(observed, symmetric_template)

        # C8BC asymmetric model is preferred if its matched filter gives
        # a higher statistic than the symmetric template
        if stat_c8bc > stat_symmetric:
            correct_identifications += 1

    return correct_identifications / n_realizations


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_results(amplitudes, probabilities, uncertainties,
                  false_positive_rate, asymmetry_prob,
                  save_path="c8bc_northern_cold_spot_detection.png"):
    """
    Generate the supplementary figure for Code S2.
    Three panels: detection probability curve, comparison with null,
    and structural asymmetry discrimination probability.
    """
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.35)

    # Panel 1 — Detection probability vs amplitude
    ax1 = fig.add_subplot(gs[0])
    ax1.fill_between(amplitudes,
                      probabilities - uncertainties,
                      probabilities + uncertainties,
                      alpha=0.3, color='steelblue', label='1σ uncertainty')
    ax1.plot(amplitudes, probabilities,
             color='steelblue', linewidth=2, label='C8BC signal')
    ax1.axhline(y=0.95, color='darkgreen', linestyle='--',
                linewidth=1.5, label='95% detection threshold')
    ax1.axhline(y=0.50, color='gray', linestyle=':', linewidth=1)
    ax1.axvline(x=SIGNAL_AMPLITUDE_CENTRAL_UK, color='orange',
                linestyle='--', linewidth=1.5, label=f'Central estimate '
                f'({SIGNAL_AMPLITUDE_CENTRAL_UK} μK)')
    ax1.set_xlabel('Northern Cold Spot Amplitude (μK)', fontsize=11)
    ax1.set_ylabel('Detection Probability', fontsize=11)
    ax1.set_title('CMB-S4 Detection Probability\nvs Signal Amplitude',
                   fontsize=11)
    ax1.legend(fontsize=8)
    ax1.set_ylim(0, 1.05)
    ax1.grid(True, alpha=0.3)

    # Panel 2 — C8BC signal vs null false positive rate
    ax2 = fig.add_subplot(gs[1])
    central_idx = np.argmin(np.abs(amplitudes - SIGNAL_AMPLITUDE_CENTRAL_UK))
    central_prob = probabilities[central_idx]
    central_unc = uncertainties[central_idx]

    categories = ['C8BC\n(central estimate)', 'Null hypothesis\n(false positive)']
    values = [central_prob, false_positive_rate]
    errors = [central_unc, np.sqrt(false_positive_rate *
                                    (1 - false_positive_rate) / N_MONTE_CARLO)]
    colors = ['steelblue', 'tomato']

    bars = ax2.bar(categories, values, color=colors, alpha=0.8,
                   yerr=errors, capsize=6)
    ax2.set_ylabel('Probability', fontsize=11)
    ax2.set_title(f'Signal vs Null at {DETECTION_THRESHOLD_SIGMA}σ\n'
                   f'Detection Threshold', fontsize=11)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width()/2., val + 0.02,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=10,
                 fontweight='bold')

    # Panel 3 — Structural asymmetry discrimination
    ax3 = fig.add_subplot(gs[2])
    asymmetry_unc = np.sqrt(
        asymmetry_prob * (1 - asymmetry_prob) / N_MONTE_CARLO
    )
    ax3.bar(['C8BC asymmetric\nmodel preferred',
              'Symmetric model\npreferred'],
             [asymmetry_prob, 1 - asymmetry_prob],
             color=['steelblue', 'tomato'], alpha=0.8,
             yerr=[asymmetry_unc, asymmetry_unc], capsize=6)
    ax3.axhline(y=0.5, color='gray', linestyle='--',
                linewidth=1.5, label='Random chance')
    ax3.set_ylabel('Probability', fontsize=11)
    ax3.set_title('Structural Asymmetry\nDiscrimination Probability',
                   fontsize=11)
    ax3.set_ylim(0, 1.05)
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.text(0, asymmetry_prob + 0.02, f'{asymmetry_prob:.3f}',
              ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.suptitle('C8BC Northern Cold Spot: CMB-S4 Detection Analysis\n'
                  f'(N = {N_MONTE_CARLO:,} Monte Carlo realizations, '
                  f'noise = {NOISE_UK_PER_SQRDEG} μK/deg², '
                  f'threshold = {DETECTION_THRESHOLD_SIGMA}σ)',
                  fontsize=10, y=1.02)

    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {save_path}")
    return fig


# ── Main execution ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("C8BC Supplementary Code S2")
    print("Northern Cold Spot Detection Probability Analysis")
    print("=" * 65)
    print()

    # Step 1 — Amplitude scan
    print("Step 1: Scanning detection probability over C8BC amplitude range")
    amplitudes, probabilities, uncertainties = amplitude_scan(n_points=20)

    # Step 2 — Null hypothesis false positive rate
    print("\nStep 2: Computing null hypothesis false positive rate")
    false_positive_rate = compute_null_false_positive_rate()
    print(f"  False positive rate at {DETECTION_THRESHOLD_SIGMA}σ: "
          f"{false_positive_rate:.4f}")

    # Step 3 — Structural asymmetry discrimination
    print("\nStep 3: Computing structural asymmetry discrimination probability")
    asymmetry_prob = structural_asymmetry_discriminant()
    print(f"  P(C8BC asymmetric model preferred over symmetric): "
          f"{asymmetry_prob:.3f}")

    # Step 4 — Summary statistics
    central_idx = np.argmin(
        np.abs(amplitudes - SIGNAL_AMPLITUDE_CENTRAL_UK)
    )
    central_prob = probabilities[central_idx]
    min_detectable = amplitudes[probabilities >= 0.95][0] \
        if any(probabilities >= 0.95) else ">140 μK"

    print("\n" + "=" * 65)
    print("SUMMARY RESULTS")
    print("=" * 65)
    print(f"C8BC central estimate ({SIGNAL_AMPLITUDE_CENTRAL_UK} μK):")
    print(f"  Detection probability: {central_prob:.3f} "
          f"+/- {uncertainties[central_idx]:.3f}")
    print(f"  False positive rate:   {false_positive_rate:.4f}")
    print(f"  Signal-to-noise ratio: "
          f"{central_prob/max(false_positive_rate, 1e-6):.1f}x over null")
    print(f"Minimum detectable amplitude at 95% probability: "
          f"{min_detectable} μK")
    print(f"Structural asymmetry discrimination: {asymmetry_prob:.3f}")
    print(f"  (0.5 = random chance, 1.0 = perfect discrimination)")
    print("=" * 65)

    # Step 5 — Generate figure
    print("\nStep 5: Generating supplementary figure")
    plot_results(amplitudes, probabilities, uncertainties,
                  false_positive_rate, asymmetry_prob)

    print("\nCode S2 complete.")
    print("Cite as: Hurst (2026), C8BC Supplementary Code S2,")
    print("GitHub: https://github.com/Aaronius0000/C8BC-MonteCarlo")
