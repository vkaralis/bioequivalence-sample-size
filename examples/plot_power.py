"""Create a power-versus-sample-size graph (requires the 'plot' extra)."""

import matplotlib.pyplot as plt

from sample_size import power_crossover_log

sample_sizes = list(range(12, 38, 2))
powers = power_crossover_log(sample_sizes, gmr=0.90, cv_within=0.20)

plt.plot(sample_sizes, [100 * value for value in powers], "o-")
plt.xlabel("Total sample size")
plt.ylabel("Approximate power (%)")
plt.ylim(0, 100)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("power_curve.png", dpi=160)

