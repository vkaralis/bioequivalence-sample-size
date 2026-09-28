"""Command-line interface for the sample-size package."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from .core import (
    crossover_additive,
    crossover_log,
    crossover_log_iterative,
    parallel_additive,
    parallel_log,
    power_crossover_log,
)


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--power", type=float, default=0.80)
    parser.add_argument("--no-correction", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Equivalence-study sample-size calculator")
    sub = parser.add_subparsers(dest="command", required=True)

    pa = sub.add_parser("parallel-additive")
    pa.add_argument("--mean-test", type=float, required=True)
    pa.add_argument("--mean-reference", type=float, required=True)
    pa.add_argument("--margin", type=float, required=True)
    pa.add_argument("--sd", type=float, required=True)
    pa.add_argument("--ratio", type=float, default=1.0)
    _common(pa)

    pl = sub.add_parser("parallel-log")
    pl.add_argument("--gmr", type=float, required=True)
    pl.add_argument("--sd-log", type=float, required=True)
    pl.add_argument("--lower", type=float, default=0.80)
    pl.add_argument("--upper", type=float, default=1.25)
    pl.add_argument("--ratio", type=float, default=1.0)
    _common(pl)

    ca = sub.add_parser("crossover-additive")
    ca.add_argument("--mean-test", type=float, required=True)
    ca.add_argument("--mean-reference", type=float, required=True)
    ca.add_argument("--margin", type=float, required=True)
    ca.add_argument("--sd-within", type=float, required=True)
    _common(ca)

    for name in ("crossover-log", "crossover-iterative"):
        cl = sub.add_parser(name)
        cl.add_argument("--gmr", type=float, required=True)
        variability = cl.add_mutually_exclusive_group(required=True)
        variability.add_argument("--cv-within", type=float)
        variability.add_argument("--sd-log", type=float)
        cl.add_argument("--lower", type=float, default=0.80)
        cl.add_argument("--upper", type=float, default=1.25)
        if name == "crossover-iterative":
            cl.add_argument("--max-sample-size", type=int, default=10_000)
        _common(cl)

    pw = sub.add_parser("power")
    pw.add_argument("--n", type=int, nargs="+", required=True)
    pw.add_argument("--gmr", type=float, required=True)
    variability = pw.add_mutually_exclusive_group(required=True)
    variability.add_argument("--cv-within", type=float)
    variability.add_argument("--sd-log", type=float)
    pw.add_argument("--lower", type=float, default=0.80)
    pw.add_argument("--upper", type=float, default=1.25)
    pw.add_argument("--alpha", type=float, default=0.05)
    pw.add_argument("--no-correction", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    correction = not args.no_correction
    if args.command == "parallel-additive":
        result = parallel_additive(args.mean_test, args.mean_reference, args.margin, args.sd,
                                   alpha=args.alpha, power=args.power,
                                   allocation_ratio=args.ratio, correction=correction)
    elif args.command == "parallel-log":
        result = parallel_log(args.gmr, args.sd_log, lower=args.lower, upper=args.upper,
                              alpha=args.alpha, power=args.power,
                              allocation_ratio=args.ratio, correction=correction)
    elif args.command == "crossover-additive":
        result = crossover_additive(args.mean_test, args.mean_reference, args.margin,
                                    args.sd_within, alpha=args.alpha, power=args.power,
                                    correction=correction)
    elif args.command in {"crossover-log", "crossover-iterative"}:
        function = crossover_log_iterative if args.command.endswith("iterative") else crossover_log
        options = {
            "cv_within": args.cv_within,
            "sd_log": args.sd_log,
            "lower": args.lower,
            "upper": args.upper,
            "alpha": args.alpha,
            "power": args.power,
            "correction": correction,
        }
        if args.command == "crossover-iterative":
            options["max_sample_size"] = args.max_sample_size
        result = function(args.gmr, **options)
    else:
        powers = power_crossover_log(args.n, args.gmr, cv_within=args.cv_within,
                                     sd_log=args.sd_log, lower=args.lower, upper=args.upper,
                                     alpha=args.alpha, correction=correction)
        print(json.dumps({"sample_sizes": args.n, "power": powers}, indent=2))
        return
    print(json.dumps(asdict(result), indent=2))


if __name__ == "__main__":
    main()
