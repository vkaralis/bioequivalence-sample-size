"""Sample-size formulae for equivalence studies.

The functions implement normal approximations for equivalence studies, the
Guenther correction described by Julious (2004), and a discrete power search
for a balanced 2x2 crossover design.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import ceil, exp, isfinite, log, sqrt

from scipy.stats import norm, t


@dataclass(frozen=True)
class SampleSizeResult:
    """A sample-size result before and after design-compatible rounding."""

    raw: float
    total: int
    method: str
    n_test: int | None = None
    n_reference: int | None = None
    iterations: int = 0


def _positive(name: str, value: float) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite number greater than zero")


def _probabilities(alpha: float, power: float) -> None:
    if not 0 < alpha < 0.5:
        raise ValueError("alpha must be between 0 and 0.5 (one-sided TOST alpha)")
    if not 0.5 < power < 1:
        raise ValueError("power must be between 0.5 and 1")


def _validate_limits(lower: float, upper: float) -> None:
    """Validate positive, ordered multiplicative equivalence limits."""
    if not isfinite(lower) or not isfinite(upper) or lower <= 0 or upper <= 0:
        raise ValueError("equivalence limits must be finite numbers greater than zero")
    if lower >= upper:
        raise ValueError("lower equivalence limit must be smaller than upper")


def _margin(estimate: float, lower: float, upper: float) -> float:
    if not lower < estimate < upper:
        raise ValueError("the expected effect must lie strictly inside the equivalence limits")
    return min(estimate - lower, upper - estimate)


def _round_even(value: float) -> int:
    result = ceil(value)
    return result if result % 2 == 0 else result + 1


def _z_sum(alpha: float, power: float, centered: bool = False) -> tuple[float, float]:
    za = float(norm.ppf(1 - alpha))
    # At the exact centre both equivalence boundaries contribute symmetrically,
    # matching the convention in the supplied iterative MATLAB implementation.
    beta_tail = (1 - power) / 2 if centered else 1 - power
    return za, float(norm.ppf(1 - beta_tail))


def cv_to_log_sd(cv: float) -> float:
    """Convert a within-subject coefficient of variation to log-scale SD."""
    _positive("cv", cv)
    return sqrt(log(1 + cv**2))


def _parallel(
    estimate: float,
    lower: float,
    upper: float,
    sd: float,
    alpha: float,
    power: float,
    allocation_ratio: float,
    correction: bool,
    method: str,
) -> SampleSizeResult:
    _positive("sd", sd)
    _positive("allocation_ratio", allocation_ratio)
    _probabilities(alpha, power)
    margin = _margin(estimate, lower, upper)
    centered = abs(estimate - (lower + upper) / 2) < 1e-12
    za, zb = _z_sum(alpha, power, centered)
    n_reference_raw = ((allocation_ratio + 1) / allocation_ratio) * sd**2 * (
        za + zb
    ) ** 2 / margin**2
    if correction:
        # Guenther correction for the per-reference-group sample size; see
        # Julious (2004), equation 2.8. For equal allocation this adds z²/4
        # to each arm (z²/2 to the total sample size).
        n_reference_raw += za**2 / 4
    n_reference = ceil(n_reference_raw)
    n_test = ceil(allocation_ratio * n_reference)
    return SampleSizeResult(
        raw=n_reference_raw * (1 + allocation_ratio),
        total=n_reference + n_test,
        n_test=n_test,
        n_reference=n_reference,
        method=method + (" with correction" if correction else ""),
    )


def parallel_additive(
    mean_test: float,
    mean_reference: float,
    equivalence_margin: float,
    sd: float,
    *,
    alpha: float = 0.05,
    power: float = 0.80,
    allocation_ratio: float = 1.0,
    correction: bool = True,
) -> SampleSizeResult:
    """Sample size for a parallel, additive-model equivalence study."""
    _positive("equivalence_margin", equivalence_margin)
    difference = mean_test - mean_reference
    return _parallel(
        difference,
        -equivalence_margin,
        equivalence_margin,
        sd,
        alpha,
        power,
        allocation_ratio,
        correction,
        "parallel additive",
    )


def parallel_log(
    gmr: float,
    sd_log: float,
    *,
    lower: float = 0.80,
    upper: float = 1.25,
    alpha: float = 0.05,
    power: float = 0.80,
    allocation_ratio: float = 1.0,
    correction: bool = True,
) -> SampleSizeResult:
    """Sample size for a parallel multiplicative model using log-scale SD."""
    _positive("gmr", gmr)
    _validate_limits(lower, upper)
    return _parallel(
        log(gmr),
        log(lower),
        log(upper),
        sd_log,
        alpha,
        power,
        allocation_ratio,
        correction,
        "parallel log",
    )


def _crossover(
    estimate: float,
    lower: float,
    upper: float,
    sd: float,
    alpha: float,
    power: float,
    correction: bool,
    method: str,
) -> SampleSizeResult:
    _positive("sd", sd)
    _probabilities(alpha, power)
    margin = _margin(estimate, lower, upper)
    centered = abs(estimate - (lower + upper) / 2) < 1e-12
    za, zb = _z_sum(alpha, power, centered)
    raw = 2 * sd**2 * (za + zb) ** 2 / margin**2
    if correction:
        # Guenther correction for total n in a balanced crossover; Julious
        # (2004), equation 2.14.
        raw += za**2 / 2
    total = _round_even(raw)
    return SampleSizeResult(raw, total, method + (" with correction" if correction else ""))


def crossover_additive(
    mean_test: float,
    mean_reference: float,
    equivalence_margin: float,
    sd_within: float,
    *,
    alpha: float = 0.05,
    power: float = 0.80,
    correction: bool = True,
) -> SampleSizeResult:
    """Sample size for a balanced 2x2 crossover additive model."""
    _positive("equivalence_margin", equivalence_margin)
    return _crossover(
        mean_test - mean_reference,
        -equivalence_margin,
        equivalence_margin,
        sd_within,
        alpha,
        power,
        correction,
        "2x2 crossover additive",
    )


def crossover_log(
    gmr: float,
    *,
    cv_within: float | None = None,
    sd_log: float | None = None,
    lower: float = 0.80,
    upper: float = 1.25,
    alpha: float = 0.05,
    power: float = 0.80,
    correction: bool = True,
) -> SampleSizeResult:
    """Sample size for a balanced 2x2 crossover multiplicative model."""
    sd = _resolve_log_sd(cv_within, sd_log)
    _positive("gmr", gmr)
    _validate_limits(lower, upper)
    return _crossover(
        log(gmr), log(lower), log(upper), sd, alpha, power, correction, "2x2 crossover log"
    )


def _resolve_log_sd(cv_within: float | None, sd_log: float | None) -> float:
    if (cv_within is None) == (sd_log is None):
        raise ValueError("provide exactly one of cv_within or sd_log")
    if cv_within is not None:
        return cv_to_log_sd(cv_within)
    assert sd_log is not None
    _positive("sd_log", sd_log)
    return sd_log


def crossover_log_iterative(
    gmr: float,
    *,
    cv_within: float | None = None,
    sd_log: float | None = None,
    lower: float = 0.80,
    upper: float = 1.25,
    alpha: float = 0.05,
    power: float = 0.80,
    correction: bool = True,
    max_sample_size: int = 10_000,
) -> SampleSizeResult:
    """Find the smallest even sample size attaining the approximate power.

    The search evaluates every admissible balanced total sample size (4, 6,
    8, ...) and therefore cannot oscillate like a fixed-point iteration.
    """
    _resolve_log_sd(cv_within, sd_log)
    _positive("gmr", gmr)
    _probabilities(alpha, power)
    _validate_limits(lower, upper)
    _margin(log(gmr), log(lower), log(upper))
    if isinstance(max_sample_size, bool) or not isinstance(max_sample_size, int):
        raise ValueError(  # noqa: TRY004 - public validation consistently uses ValueError
            "max_sample_size must be an even integer of at least 4"
        )
    if max_sample_size < 4 or max_sample_size % 2:
        raise ValueError("max_sample_size must be an even integer of at least 4")

    # At n=4, t_(0.95, 2)^2 / 2 exceeds n, so the corrected approximation
    # has no positive effective sample size. Begin at the first evaluable n.
    first_sample_size = 6 if correction else 4
    for iteration, sample_size in enumerate(
        range(first_sample_size, max_sample_size + 1, 2), start=1
    ):
        achieved = power_crossover_log(
            sample_size,
            gmr,
            cv_within=cv_within,
            sd_log=sd_log,
            lower=lower,
            upper=upper,
            alpha=alpha,
            correction=correction,
        )
        if achieved >= power:
            return SampleSizeResult(
                float(sample_size),
                sample_size,
                "2x2 crossover log power search"
                + (" with Guenther correction" if correction else ""),
                iterations=iteration,
            )
    raise RuntimeError(f"target power was not reached by total sample size {max_sample_size}")


def power_crossover_log(
    sample_sizes: int | Iterable[int],
    gmr: float,
    *,
    cv_within: float | None = None,
    sd_log: float | None = None,
    lower: float = 0.80,
    upper: float = 1.25,
    alpha: float = 0.05,
    correction: bool = True,
) -> float | list[float]:
    """Approximate power for balanced 2x2 crossover sample sizes.

    Returns proportions in [0, 1]. This is the central-t approximation in
    Julious (2004), equation 6.2, evaluated at both equivalence boundaries.
    With ``correction=True``, the inverse Guenther adjustment n - t²/2 is used.
    """
    sd = _resolve_log_sd(cv_within, sd_log)
    _positive("gmr", gmr)
    _validate_limits(lower, upper)
    if not 0 < alpha < 0.5:
        raise ValueError("alpha must be between 0 and 0.5")
    theta = log(gmr)
    lower_log = log(lower)
    upper_log = log(upper)
    _margin(theta, lower_log, upper_log)
    if isinstance(sample_sizes, bool):
        raise ValueError(  # noqa: TRY004 - bool is invalid as a numeric sample size
            "sample_sizes must contain even integers of at least 4"
        )
    scalar = isinstance(sample_sizes, int)
    if scalar:
        values = [sample_sizes]
    else:
        try:
            values = list(sample_sizes)
        except TypeError as error:
            raise ValueError("sample_sizes must be an integer or an iterable of integers") from error
    if not values:
        raise ValueError("sample_sizes must contain at least one sample size")
    powers: list[float] = []
    for n in values:
        if isinstance(n, bool) or not isinstance(n, int) or n < 4 or n % 2:
            raise ValueError(
                "each sample size must be an even integer of at least 4 "
                "for a balanced 2x2 crossover"
            )
        critical = float(t.ppf(1 - alpha, n - 2))
        effective_n = n - critical**2 / 2 if correction else float(n)
        if effective_n <= 0:
            raise ValueError(
                f"effective sample size is not positive for n={n}; "
                "increase n or disable the Guenther correction"
            )
        scale = sqrt(effective_n / (2 * sd**2))
        lower_score = scale * (theta - lower_log) - critical
        upper_score = scale * (upper_log - theta) - critical
        approximate_power = norm.cdf(lower_score) + norm.cdf(upper_score) - 1
        powers.append(min(1.0, max(0.0, float(approximate_power))))
    return powers[0] if scalar else powers


def cv_from_log_sd(sd_log: float) -> float:
    """Convert log-scale SD to coefficient of variation."""
    _positive("sd_log", sd_log)
    return sqrt(exp(sd_log**2) - 1)
