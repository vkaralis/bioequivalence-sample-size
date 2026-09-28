# 2×2 Bioequivalence Sample Size

Python εργαλείο για τον υπολογισμό μεγέθους δείγματος σε μελέτες
βιοϊσοδυναμίας 2×2 crossover, με συμπληρωματικές συναρτήσεις για parallel
designs και εκτίμηση ισχύος. Το project αποτελεί καθαρή, ελεγμένη μεταφορά των
MATLAB scripts του repository και υποστηρίζει:

**Topics:** `bioequivalence` · `sample-size-estimation` · `crossover` · `parallel`

- παράλληλο σχέδιο με μη μετασχηματισμένα δεδομένα (προσθετικό μοντέλο),
- παράλληλο σχέδιο με λογαριθμικά μετασχηματισμένα δεδομένα,
- ισορροπημένο 2×2 crossover με προσθετικό μοντέλο,
- ισορροπημένο 2×2 crossover με λογαριθμικό μοντέλο και `CVw` ή log-SD,
- διακριτή αναζήτηση του ελάχιστου άρτιου crossover μεγέθους με Student-t critical value,
- προσεγγιστική ισχύ για δεδομένο συνολικό μέγεθος δείγματος.

> **Σημαντικό:** Το λογισμικό προορίζεται για ερευνητική/εκπαιδευτική χρήση.
> Πριν χρησιμοποιηθεί σε πρωτόκολλο ή κανονιστική υποβολή, τα αποτελέσματα πρέπει
> να επιβεβαιώνονται από αρμόδιο βιοστατιστικό και με επικυρωμένο λογισμικό.

## Εγκατάσταση

Απαιτείται Python 3.9 ή νεότερη.

```bash
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -e .
```

Για tests και διάγραμμα ισχύος:

```bash
python -m pip install -e ".[dev,plot]"
pytest
```

## Γρήγορη χρήση από Python

```python
from sample_size import crossover_log, crossover_log_iterative

approximate = crossover_log(
    gmr=1.05,
    cv_within=0.20,
    lower=0.80,
    upper=1.25,
    alpha=0.05,       # μονόπλευρο alpha για κάθε TOST έλεγχο
    power=0.90,
)

iterative = crossover_log_iterative(
    gmr=1.05,
    cv_within=0.20,
    power=0.90,
)

print(approximate.total)
print(iterative.total, iterative.iterations)
```

Το `total` είναι ήδη στρογγυλοποιημένο σε άρτιο αριθμό για ισομερή κατανομή
στις δύο ακολουθίες crossover. Στις κλειστές προσεγγίσεις το `raw` κρατά τη μη
στρογγυλοποιημένη εκτίμηση. Στη διακριτή power search είναι ίσο με το επιλεγμένο
άθροισμα συμμετεχόντων.

### Όλες οι περιπτώσεις

```python
from sample_size import (
    parallel_additive,
    parallel_log,
    crossover_additive,
    crossover_log,
    power_crossover_log,
)

# Parallel, additive
r1 = parallel_additive(5, 4, equivalence_margin=1.25, sd=0.1)

# Parallel, log scale; ratio = n_test / n_reference
r2 = parallel_log(gmr=1.05, sd_log=0.25, allocation_ratio=1)

# 2x2 crossover, additive
r3 = crossover_additive(5, 4, equivalence_margin=1.2, sd_within=0.5)

# 2x2 crossover, log scale. Δώστε είτε CVw είτε log-SD, όχι και τα δύο.
r4 = crossover_log(gmr=1.06, cv_within=0.32)
r5 = crossover_log(gmr=1.05, sd_log=0.25)

# Approximate power as a proportion (0–1)
powers = power_crossover_log([12, 14, 16, 18], gmr=0.90, cv_within=0.20)
```

Περισσότερα παραδείγματα υπάρχουν στα [`examples/`](examples).

## Χρήση από γραμμή εντολών

Μετά την εγκατάσταση είναι διαθέσιμη η εντολή `sample-size`:

```bash
sample-size crossover-log --gmr 1.05 --cv-within 0.20 --power 0.90
sample-size crossover-iterative --gmr 1.05 --cv-within 0.20 --power 0.90
sample-size parallel-log --gmr 1.05 --sd-log 0.25 --ratio 1
sample-size power --n 12 14 16 18 --gmr 0.90 --cv-within 0.20
```

Η έξοδος είναι JSON, ώστε να μπορεί να χρησιμοποιηθεί εύκολα σε scripts.
Δείτε όλες τις επιλογές με `sample-size --help` ή
`sample-size crossover-log --help`.

## Παραδοχές και ορισμοί

- Το `alpha` είναι το μονόπλευρο επίπεδο σημαντικότητας κάθε ελέγχου TOST.
  Για το συνηθισμένο 90% CI χρησιμοποιήστε `alpha=0.05`.
- `gmr` είναι ο αναμενόμενος γεωμετρικός λόγος Test/Reference.
- `cv_within` δίνεται ως αναλογία (`0.20`, όχι `20`) και μετατρέπεται με
  `sqrt(log(1 + CVw²))`.
- `sd_log` είναι η υπολειμματική τυπική απόκλιση στη λογαριθμική κλίμακα.
- Ένα balanced 2×2 crossover απαιτεί άρτιο συνολικό `n`, τουλάχιστον 4, ώστε οι
  συμμετέχοντες να κατανέμονται ισόποσα στις ακολουθίες TR και RT. Η power
  function απορρίπτει μονά ή μικρότερα μεγέθη. Δεν περιλαμβάνεται προσαύξηση
  για dropouts.
- Στο παράλληλο σχέδιο `allocation_ratio = n_test / n_reference`.
- Η προεπιλεγμένη διόρθωση είναι η Guenther correction όπως παρουσιάζεται από
  τον Julious (2004): `z²/4` στο ανά-reference-arm μέγεθος parallel design και
  `z²/2` στο συνολικό μέγεθος balanced crossover. Δεν χρησιμοποιείται πλέον ο
  μη τεκμηριωμένος όρος `z/2` των αρχικών MATLAB scripts. Απενεργοποιείται με
  `correction=False` ή `--no-correction`.
- Στην αντίστροφη power προσέγγιση η ίδια διόρθωση εφαρμόζεται συνεπώς ως
  `effective_n = n - t²/2`. Για `n=4` αυτό μπορεί να μην είναι θετικό· τότε
  επιστρέφεται σαφές `ValueError`.
- Η `power_crossover_log()` είναι central-t/Normal approximation που αξιολογεί
  και τα δύο όρια ισοδυναμίας. Δεν είναι η ακριβής μέθοδος noncentral-t του
  Owen και μπορεί να δώσει συντηρητικότερο αποτέλεσμα.
- Η `crossover_log_iterative()` δεν κάνει πλέον fixed-point iteration και δεν
  χρησιμοποιεί `round()`. Ελέγχει διαδοχικά τα επιτρεπτά άρτια μεγέθη μέχρι να
  επιτευχθεί η ζητούμενη προσεγγιστική ισχύς ή το `max_sample_size`.

## Στατιστικές αναφορές και validation

- Julious SA. *Tutorial in Biostatistics: Sample sizes for clinical trials with
  Normal data*. Statistics in Medicine. 2004;23:1921–1986.
  [doi:10.1002/sim.1783](https://doi.org/10.1002/sim.1783).
- Το regression test `CVw=20%`, `GMR=1.00`, όρια `0.80–1.25`, power `90%`
  συγκρίνεται με τον Πίνακα 6.1 του Julious: ο πίνακας δίνει `n=19` με
  noncentral-t και το project επιστρέφει `n=20` μετά την απαραίτητη εξισορρόπηση.

## Έλεγχοι

```bash
pytest
ruff check .
```

Το GitHub Actions workflow εκτελεί αυτόματα τους ελέγχους σε Python 3.9–3.13.

## Δομή

```text
src/sample_size/       βιβλιοθήκη και CLI
tests/                 αυτοματοποιημένοι έλεγχοι
examples/              παραδείγματα και γράφημα ισχύος
*.m                    αρχικά MATLAB scripts για ιχνηλασιμότητα
```

## Άδεια

MIT — δείτε το [LICENSE](LICENSE).
