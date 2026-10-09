"""Reproduce the primary six-cohort random-effects synthesis."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, t


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "primary_cohort_estimates.csv"
DEFAULT_EXPECTED = ROOT / "data" / "reported_primary_pooled_results.csv"
DEFAULT_OUTPUT = ROOT / "output" / "reproduced_primary_pooled_results.csv"

LABELS = {
    "first_disability": "First incident ADL disability; death is a competing event",
    "death_before_disability": "Death before incident ADL disability",
    "first_exit": "First disability or death from baseline independence",
}


def reml_meta(frame: pd.DataFrame) -> dict[str, float | int]:
    """Pool log IRRs using REML and modified Hartung-Knapp inference."""
    data = frame[["log_irr_mi", "se"]].dropna()
    yi = data["log_irr_mi"].to_numpy(float)
    vi = np.square(data["se"].to_numpy(float))
    k = len(yi)
    if k < 2:
        raise ValueError("At least two cohort estimates are required")

    def objective(tau2: float) -> float:
        weights = 1.0 / (vi + tau2)
        mean = float(np.sum(weights * yi) / np.sum(weights))
        return float(
            np.sum(np.log(vi + tau2))
            + np.log(np.sum(weights))
            + np.sum(weights * np.square(yi - mean))
        )

    upper = max(1.0, float(np.var(yi, ddof=1) * 10 + np.max(vi)))
    optimization = minimize_scalar(objective, bounds=(0.0, upper), method="bounded")
    tau2 = max(0.0, float(optimization.x))
    weights = 1.0 / (vi + tau2)
    pooled = float(np.sum(weights * yi) / np.sum(weights))
    q_star = float(np.sum(weights * np.square(yi - pooled)))
    hk_se = math.sqrt(max(1.0, q_star / (k - 1)) / float(np.sum(weights)))
    critical = float(t.ppf(0.975, df=k - 1))
    lcl = pooled - critical * hk_se
    ucl = pooled + critical * hk_se
    prediction_se = math.sqrt(tau2 + hk_se**2)
    prediction_lcl = pooled - critical * prediction_se
    prediction_ucl = pooled + critical * prediction_se

    fixed_weights = 1.0 / vi
    fixed_mean = float(np.sum(fixed_weights * yi) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(yi - fixed_mean)))
    q_df = k - 1
    i2 = max(0.0, 100.0 * (q - q_df) / q) if q > 0 else 0.0

    return {
        "k": k,
        "pooled_IRR": math.exp(pooled),
        "LCL": math.exp(lcl),
        "UCL": math.exp(ucl),
        "prediction_LCL": math.exp(prediction_lcl),
        "prediction_UCL": math.exp(prediction_ucl),
        "I2_percent": i2,
        "tau2_REML": tau2,
        "p_pooled": float(2 * t.sf(abs(pooled / hk_se), df=k - 1)),
        "p_heterogeneity": float(chi2.sf(q, q_df)),
    }


def reproduce(input_path: Path) -> pd.DataFrame:
    estimates = pd.read_csv(input_path)
    required = {"cohort", "outcome", "log_irr_mi", "se"}
    missing = required.difference(estimates.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    rows = []
    for outcome in ("death_before_disability", "first_disability", "first_exit"):
        subset = estimates.loc[estimates["outcome"].eq(outcome)].copy()
        result = reml_meta(subset)
        rows.append(
            {
                "level": "6 cohorts",
                "estimand": outcome,
                "estimand_label": LABELS[outcome],
                "age_model": "restricted_cubic_spline",
                "adjustment": "Multiple imputation",
                **result,
            }
        )
    return pd.DataFrame(rows)


def verify(reproduced: pd.DataFrame, expected_path: Path) -> None:
    expected = pd.read_csv(expected_path)
    numeric = [
        "pooled_IRR",
        "LCL",
        "UCL",
        "prediction_LCL",
        "prediction_UCL",
        "I2_percent",
        "tau2_REML",
        "p_pooled",
    ]
    left = reproduced.set_index("estimand").sort_index()
    right = expected.set_index("estimand").sort_index()
    if list(left.index) != list(right.index):
        raise AssertionError("Reproduced and expected estimands differ")
    for column in numeric:
        if not np.allclose(left[column], right[column], rtol=1e-7, atol=1e-9):
            difference = (left[column] - right[column]).abs().max()
            raise AssertionError(f"Mismatch in {column}; maximum absolute difference={difference}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--expected", type=Path, default=DEFAULT_EXPECTED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    result = reproduce(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    if args.check:
        verify(result, args.expected)
        print("Primary pooled estimates reproduced and verified.")
    else:
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

