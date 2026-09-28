import math

import pytest

from sample_size import (
    crossover_additive,
    crossover_log,
    crossover_log_iterative,
    cv_from_log_sd,
    cv_to_log_sd,
    parallel_additive,
    parallel_log,
    power_crossover_log,
)


def test_cv_conversions_are_numerical_inverses() -> None:
    cv = 0.20
    sd = cv_to_log_sd(cv)
    assert sd == pytest.approx(math.sqrt(math.log(1.04)))
    assert cv_from_log_sd(sd) == pytest.approx(cv)


@pytest.mark.parametrize(
    ("correction", "raw", "total", "per_arm"),
    [(False, 3.9568366284926513, 4, 2), (True, 5.309608355540357, 6, 3)],
)
def test_parallel_additive_reference_values(
    correction: bool, raw: float, total: int, per_arm: int
) -> None:
    result = parallel_additive(5, 4, 1.25, 0.1, correction=correction)
    assert result.raw == pytest.approx(raw)
    assert result.total == total
    assert result.n_test == result.n_reference == per_arm


@pytest.mark.parametrize(
    ("correction", "raw", "total", "per_arm"),
    [(False, 50.844897530288414, 52, 26), (True, 52.197669257336116, 54, 27)],
)
def test_parallel_log_balanced_reference_values(
    correction: bool, raw: float, total: int, per_arm: int
) -> None:
    result = parallel_log(1.05, 0.25, correction=correction)
    assert result.raw == pytest.approx(raw)
    assert result.total == total
    assert result.n_test == result.n_reference == per_arm


def test_parallel_log_unequal_allocation_reference_value() -> None:
    result = parallel_log(1.05, 0.25, allocation_ratio=2)
    assert result.raw == pytest.approx(59.229667312146006)
    assert result.n_reference == 20
    assert result.n_test == 40
    assert result.total == 60


def test_crossover_correction_uses_squared_critical_value() -> None:
    uncorrected = crossover_additive(5, 4, 1.2, 0.5, correction=False)
    corrected = crossover_additive(5, 4, 1.2, 0.5, correction=True)
    z = 1.6448536269514722
    assert corrected.raw - uncorrected.raw == pytest.approx(z**2 / 2)
    assert corrected.total % 2 == 0


def test_published_julious_crossover_reference_case() -> None:
    # Julious (2004), Table 6.1: CV=20%, GMR=1.00, 80-125%, 90% power
    # gives total n=19 using the paper's noncentral-t method. A balanced design
    # requires rounding to 20, which the discrete approximation also returns.
    result = crossover_log_iterative(1.00, cv_within=0.20, power=0.90)
    assert result.total == 20


def test_low_cv_case_does_not_oscillate_or_fail_to_converge() -> None:
    result = crossover_log_iterative(0.95, cv_within=0.05, power=0.90)
    assert result.total == 6
    assert result.iterations == 1
    assert power_crossover_log(6, 0.95, cv_within=0.05) >= 0.90


def test_power_search_returns_smallest_qualifying_even_size() -> None:
    result = crossover_log_iterative(1.05, cv_within=0.20, power=0.90)
    assert result.total == 26
    assert power_crossover_log(24, 1.05, cv_within=0.20) < 0.90
    assert power_crossover_log(26, 1.05, cv_within=0.20) >= 0.90


@pytest.mark.parametrize(
    ("correction", "expected"),
    [
        (False, [0.34884417472439777, 0.6342048165885532, 0.7973889953027808]),
        (True, [0.2982467750290576, 0.6097430988568884, 0.7828111484880327]),
    ],
)
def test_power_reference_values(correction: bool, expected: list[float]) -> None:
    values = power_crossover_log([12, 24, 36], 0.90, cv_within=0.20, correction=correction)
    assert values == pytest.approx(expected)
    assert values == sorted(values)


@pytest.mark.parametrize("sample_sizes", [2, 3, 5, 7, True, [4, 5], []])
def test_power_rejects_invalid_balanced_sample_sizes(sample_sizes: object) -> None:
    with pytest.raises(ValueError, match="sample size|sample_sizes|even integer"):
        power_crossover_log(sample_sizes, 1.0, cv_within=0.20)  # type: ignore[arg-type]


def test_corrected_power_rejects_nonpositive_effective_n() -> None:
    with pytest.raises(ValueError, match="effective sample size is not positive"):
        power_crossover_log(4, 1.0, cv_within=0.20, correction=True)


@pytest.mark.parametrize(
    ("lower", "upper"),
    [(0, 1.25), (-0.8, 1.25), (1.25, 1.25), (1.30, 1.25)],
)
@pytest.mark.parametrize("function_name", ["parallel", "crossover", "iterative", "power"])
def test_invalid_limits_raise_specific_value_error(
    lower: float, upper: float, function_name: str
) -> None:
    with pytest.raises(ValueError, match="equivalence limit"):
        if function_name == "parallel":
            parallel_log(1.0, 0.20, lower=lower, upper=upper)
        elif function_name == "crossover":
            crossover_log(1.0, cv_within=0.20, lower=lower, upper=upper)
        elif function_name == "iterative":
            crossover_log_iterative(1.0, cv_within=0.20, lower=lower, upper=upper)
        else:
            power_crossover_log(12, 1.0, cv_within=0.20, lower=lower, upper=upper)


@pytest.mark.parametrize("gmr", [0.80, 1.25, 1.30])
@pytest.mark.parametrize("function_name", ["parallel", "crossover", "iterative", "power"])
def test_gmr_must_be_strictly_inside_limits(gmr: float, function_name: str) -> None:
    with pytest.raises(ValueError, match="strictly inside"):
        if function_name == "parallel":
            parallel_log(gmr, 0.20)
        elif function_name == "crossover":
            crossover_log(gmr, cv_within=0.20)
        elif function_name == "iterative":
            crossover_log_iterative(gmr, cv_within=0.20)
        else:
            power_crossover_log(12, gmr, cv_within=0.20)


def test_exactly_one_variability_measure_is_required() -> None:
    with pytest.raises(ValueError):
        crossover_log(1.0)
    with pytest.raises(ValueError):
        crossover_log(1.0, cv_within=0.20, sd_log=0.20)


def test_power_search_has_clear_upper_bound_failure() -> None:
    with pytest.raises(RuntimeError, match="target power was not reached"):
        crossover_log_iterative(1.20, cv_within=0.50, power=0.90, max_sample_size=20)
