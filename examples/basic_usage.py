"""Minimal examples for every supported design."""

from sample_size import (
    crossover_additive,
    crossover_log,
    crossover_log_iterative,
    parallel_additive,
    parallel_log,
    power_crossover_log,
)

print(parallel_additive(5, 4, equivalence_margin=1.25, sd=0.1))
print(parallel_log(gmr=1.05, sd_log=0.25))
print(crossover_additive(5, 4, equivalence_margin=1.2, sd_within=0.5))
print(crossover_log(gmr=1.06, cv_within=0.32, power=0.90))
print(crossover_log_iterative(gmr=1.05, cv_within=0.20, power=0.90))
print(power_crossover_log(range(12, 38, 2), gmr=0.90, cv_within=0.20))

