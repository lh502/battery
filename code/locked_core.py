#!/usr/bin/env python3
"""Route-A terminal-voltage reconstruction and endpoint-consistency audit.

The retained empirical model is

    V_hat(t) = U_dev(z_hat(t)),

where U_dev is a development-derived *loaded* terminal-voltage--SOC relation.
It is not interpreted as an OCV curve, and no R0/R1/C1 or polarization state is
identified. Because U_dev(0)=3.0 V is imposed by construction, the predicted
voltage crossing and SOC-depletion endpoint are mathematically coupled. The
endpoint-time result is therefore reported as capacity-derived endpoint timing,
not as independent validation of a voltage-triggered cutoff mechanism.

Chronological protocol
----------------------
- Cycles 1--40: initial model fitting.
- Cycles 41--60: validation audit using only the Cycles 1--40 model.
- Cycles 1--60: frozen refit after validation.
- Cycles 61--80: later-cycle robustness screening.
- Cycles 81--101: post-screening final confirmation without further changes.

Primary voltage errors are evaluated only on the common observed/predicted
support ending at the earlier endpoint. Full measured-support errors are also
saved as a diagnostic. Outputs include endpoint counts (including not reached),
a capacity-mapping diagnostic, moving-block bootstrap intervals, maximum-error
diagnostics, revised figures, and a machine-readable audit JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from sklearn.isotonic import IsotonicRegression

V_CUT = 3.0
SOC_MIN = 0.0
FIT_END = 40
VALIDATION_END = 60
ROBUSTNESS_END = 80
TIE_TOLERANCE_S = 1.0
SOC_BINS = [0.0, 0.30, 0.70, 1.0000001]
SOC_LABELS = ["0.00-0.30", "0.30-0.70", "0.70-1.00"]
BOOTSTRAP_SAMPLES = 20_000
BOOTSTRAP_SEED = 20260806
MOVING_BLOCK_LENGTH = 5


@dataclass(frozen=True)
class Profile:
    cycle: int
    time_s: np.ndarray
    current_a: np.ndarray
    voltage_v: np.ndarray
    measured_cutoff_s: float
    capacity_to_cutoff_ah: float
    n_points: int


@dataclass(frozen=True)
class EmpiricalRelation:
    soc_grid: np.ndarray
    voltage_grid: np.ndarray
    source: str

    def __call__(self, soc: np.ndarray | float) -> np.ndarray:
        z = np.asarray(soc, dtype=float)
        return np.interp(
            np.clip(z, self.soc_grid[0], self.soc_grid[-1]),
            self.soc_grid,
            self.voltage_grid,
        )


class ExperimentError(RuntimeError):
    """Raised when the data cannot support the prespecified audit."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        x = float(value)
        return x if math.isfinite(x) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def first_crossing_time(time_s: np.ndarray, values: np.ndarray, threshold: float) -> float:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    valid = np.isfinite(t) & np.isfinite(y)
    t = t[valid]
    y = y[valid]
    if len(t) == 0:
        return math.nan
    if y[0] <= threshold:
        return float(t[0])
    idx = np.where((y[:-1] > threshold) & (y[1:] <= threshold))[0]
    if len(idx) == 0:
        return math.nan
    i = int(idx[0])
    denominator = y[i] - y[i + 1]
    if abs(denominator) < 1e-15:
        return float(t[i + 1])
    fraction = (y[i] - threshold) / denominator
    return float(t[i] + fraction * (t[i + 1] - t[i]))


def interp_at(time_s: np.ndarray, values: np.ndarray, query_s: float) -> float:
    return float(np.interp(query_s, np.asarray(time_s, float), np.asarray(values, float)))


def capacity_integral_ah(time_s: np.ndarray, current_a: np.ndarray) -> float:
    if len(time_s) < 2:
        return 0.0
    return float(np.trapezoid(current_a, time_s) / 3600.0)


def arrays_to_cutoff(profile: Profile) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return measured arrays through the exact interpolated 3.0-V crossing."""
    keep = profile.time_s < profile.measured_cutoff_s - 1e-9
    t = profile.time_s[keep]
    i = profile.current_a[keep]
    v = profile.voltage_v[keep]
    t = np.append(t, profile.measured_cutoff_s)
    i = np.append(i, interp_at(profile.time_s, profile.current_a, profile.measured_cutoff_s))
    v = np.append(v, V_CUT)
    return t, i, v


def arrays_to_time(profile: Profile, end_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return measured arrays through an exact interpolated end time."""
    if not math.isfinite(end_s):
        raise ExperimentError("A finite evaluation end time is required.")
    end_s = float(np.clip(end_s, profile.time_s[0], profile.time_s[-1]))
    keep = profile.time_s < end_s - 1e-9
    t = profile.time_s[keep]
    i = profile.current_a[keep]
    v = profile.voltage_v[keep]
    t = np.append(t, end_s)
    i = np.append(i, interp_at(profile.time_s, profile.current_a, end_s))
    v = np.append(v, interp_at(profile.time_s, profile.voltage_v, end_s))
    return t, i, v


def read_profiles(raw_file: Path) -> tuple[dict[int, Profile], dict[str, Any]]:
    raw = pd.read_csv(raw_file, sep="\t")
    raw.columns = [str(column).strip() for column in raw.columns]
    required = ["Time", "mV", "mA", "Temperature", "Discharge count", "Capacity"]
    missing = [column for column in required if column not in raw.columns]
    if missing:
        raise ExperimentError(f"Missing required columns: {missing}")

    for column in required:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    rows_before = len(raw)
    raw = raw.dropna(subset=required).copy()

    discharge = raw[(raw["Discharge count"] > 0) & (raw["mA"] < -100.0)].copy()
    profiles: dict[int, Profile] = {}
    rejected: list[dict[str, Any]] = []

    for cycle_value, group in discharge.groupby("Discharge count", sort=True):
        cycle = int(cycle_value)
        group = group.sort_values("Time").drop_duplicates("Time")
        time_s = (group["Time"].to_numpy(float) - float(group["Time"].iloc[0])) * 60.0
        current_a = -group["mA"].to_numpy(float) / 1000.0
        voltage_v = group["mV"].to_numpy(float) / 1000.0
        duration_h = float(time_s[-1] / 3600.0) if len(time_s) else 0.0
        capacity_pct = float(group["Capacity"].max()) if len(group) else 0.0

        reasons: list[str] = []
        if len(group) < 80:
            reasons.append("fewer than 80 discharge samples")
        if duration_h < 1.5:
            reasons.append("duration below 1.5 h")
        if capacity_pct < 80.0:
            reasons.append("recorded discharged capacity below 80%")
        if len(time_s) < 2 or np.any(np.diff(time_s) <= 0.0):
            reasons.append("non-increasing or insufficient time stamps")
        cutoff_s = first_crossing_time(time_s, voltage_v, V_CUT)
        if not math.isfinite(cutoff_s):
            reasons.append("no measured 3.0-V first crossing")

        if reasons:
            rejected.append({"cycle": cycle, "reasons": reasons})
            continue

        keep = time_s < cutoff_s - 1e-9
        t_cut = np.append(time_s[keep], cutoff_s)
        i_cut = np.append(current_a[keep], interp_at(time_s, current_a, cutoff_s))
        capacity_ah = capacity_integral_ah(t_cut, i_cut)
        if capacity_ah <= 0.0:
            rejected.append({"cycle": cycle, "reasons": ["non-positive capacity to 3.0 V"]})
            continue

        profiles[cycle] = Profile(
            cycle=cycle,
            time_s=time_s,
            current_a=current_a,
            voltage_v=voltage_v,
            measured_cutoff_s=cutoff_s,
            capacity_to_cutoff_ah=capacity_ah,
            n_points=len(group),
        )

    if not profiles:
        raise ExperimentError("No usable cycles remain after quality control.")

    audit = {
        "raw_rows_before_required_dropna": int(rows_before),
        "raw_rows_after_required_dropna": int(len(raw)),
        "candidate_discharge_cycles": int(discharge["Discharge count"].nunique()),
        "accepted_cycles": int(len(profiles)),
        "accepted_cycle_ids": sorted(profiles),
        "rejected_cycles": rejected,
        "raw_time_unit": "minutes",
        "internal_time_unit": "seconds",
        "measured_endpoint": "linearly interpolated first crossing of measured voltage <= 3.0 V",
    }
    return profiles, audit


def fit_saturating_capacity(profiles: Iterable[Profile]) -> dict[str, float]:
    selected = sorted(list(profiles), key=lambda profile: profile.cycle)
    cycles = np.array([profile.cycle for profile in selected], dtype=float)
    capacity_mah = np.array(
        [profile.capacity_to_cutoff_ah * 1000.0 for profile in selected], dtype=float
    )
    if len(cycles) < 10:
        raise ExperimentError("Too few cycles for capacity fitting.")

    best: tuple[float, float, float, float] | None = None
    for tau in np.geomspace(3.0, 500.0, 1500):
        design = np.column_stack([np.ones_like(cycles), np.exp(-cycles / tau)])
        c_inf, amplitude = np.linalg.lstsq(design, capacity_mah, rcond=None)[0]
        if amplitude < 0.0 or not 700.0 <= c_inf <= 1000.0:
            continue
        residual = design @ np.array([c_inf, amplitude]) - capacity_mah
        sse = float(np.sum(residual**2))
        if best is None or sse < best[0]:
            best = (sse, float(c_inf), float(amplitude), float(tau))
    if best is None:
        raise ExperimentError("No feasible saturating capacity model was identified.")
    _, c_inf, amplitude, tau = best
    return {"C_inf_mah": c_inf, "amplitude_mah": amplitude, "tau_cycles": tau}


def predict_capacity_ah(cycle: int | np.ndarray, parameters: dict[str, float]) -> np.ndarray:
    n = np.asarray(cycle, dtype=float)
    capacity_mah = parameters["C_inf_mah"] + parameters["amplitude_mah"] * np.exp(
        -n / parameters["tau_cycles"]
    )
    return capacity_mah / 1000.0


def build_empirical_relation(profiles: Iterable[Profile]) -> EmpiricalRelation:
    soc_chunks: list[np.ndarray] = []
    voltage_chunks: list[np.ndarray] = []

    for profile in profiles:
        time_s, current_a, voltage_v = arrays_to_cutoff(profile)
        discharged_as = cumulative_trapezoid(current_a, time_s, initial=0.0)
        soc = 1.0 - discharged_as / (3600.0 * profile.capacity_to_cutoff_ah)
        soc[-1] = 0.0
        valid = (
            np.isfinite(soc)
            & np.isfinite(voltage_v)
            & (soc >= -1e-6)
            & (soc <= 1.0 + 1e-6)
        )
        soc_chunks.append(np.clip(soc[valid], 0.0, 1.0))
        voltage_chunks.append(voltage_v[valid])

    soc = np.concatenate(soc_chunks)
    voltage = np.concatenate(voltage_chunks)
    order = np.argsort(soc)
    soc = soc[order]
    voltage = voltage[order]

    isotonic = IsotonicRegression(increasing=True, out_of_bounds="clip")
    fitted = isotonic.fit_transform(soc, voltage)
    frame = (
        pd.DataFrame({"soc": soc, "voltage_v": fitted})
        .groupby("soc", as_index=False)
        .mean()
    )
    grid = np.linspace(0.0, 1.0, 1001)
    values = np.interp(grid, frame["soc"], frame["voltage_v"])
    values[0] = V_CUT
    values = np.maximum.accumulate(values)

    return EmpiricalRelation(
        soc_grid=grid,
        voltage_grid=values,
        source=(
            "development-derived loaded terminal-voltage--SOC relation; "
            "all construction trajectories truncated at measured 3.0-V first passage"
        ),
    )


def classify_event(t_voltage: float, t_soc: float) -> tuple[float, str]:
    voltage_finite = math.isfinite(t_voltage)
    soc_finite = math.isfinite(t_soc)
    if voltage_finite and soc_finite:
        if abs(t_voltage - t_soc) <= TIE_TOLERANCE_S:
            return min(t_voltage, t_soc), "tie"
        if t_voltage < t_soc:
            return t_voltage, "voltage-first"
        return t_soc, "SOC-first"
    if voltage_finite:
        return t_voltage, "voltage-first"
    if soc_finite:
        return t_soc, "SOC-first"
    return math.nan, "not-reached"


def error_metrics(error: np.ndarray) -> dict[str, float]:
    values = np.asarray(error, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {
            "MAE": math.nan,
            "RMSE": math.nan,
            "bias": math.nan,
            "max_abs": math.nan,
        }
    return {
        "MAE": float(np.mean(np.abs(values))),
        "RMSE": float(np.sqrt(np.mean(values**2))),
        "bias": float(np.mean(values)),
        "max_abs": float(np.max(np.abs(values))),
    }


def evaluate_window(
    window_name: str,
    profiles: list[Profile],
    capacity_parameters: dict[str, float],
    relation: EmpiricalRelation,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluate one chronological window with explicit support accounting.

    Primary voltage metrics use the common support ending at the earlier of the
    measured 3.0-V endpoint and the predicted capacity-derived endpoint. Full
    measured-support metrics are retained only as a diagnostic because the
    empirical model is clipped at U_dev(0)=3.0 V after its predicted endpoint.
    """
    cycle_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    common_support_errors: list[np.ndarray] = []
    full_support_errors: list[np.ndarray] = []
    observed_capacity_diagnostic_errors: list[np.ndarray] = []

    for profile in profiles:
        predicted_capacity_ah = float(predict_capacity_ah(profile.cycle, capacity_parameters))

        # Dense 1-s grid for endpoint timing.
        dense_time_s = np.arange(0.0, profile.time_s[-1], 1.0, dtype=float)
        if len(dense_time_s) == 0 or dense_time_s[-1] < profile.time_s[-1]:
            dense_time_s = np.append(dense_time_s, profile.time_s[-1])
        dense_current_a = np.interp(dense_time_s, profile.time_s, profile.current_a)
        discharged_as = cumulative_trapezoid(dense_current_a, dense_time_s, initial=0.0)
        dense_soc = 1.0 - discharged_as / (3600.0 * predicted_capacity_ah)
        dense_voltage = relation(dense_soc)
        t_voltage = first_crossing_time(dense_time_s, dense_voltage, V_CUT)
        t_soc = first_crossing_time(dense_time_s, dense_soc, SOC_MIN)
        predicted_endpoint_s, event_type = classify_event(t_voltage, t_soc)

        # Full observed trajectory through the exact measured 3.0-V crossing.
        full_time_s, full_current_a, measured_voltage_v = arrays_to_cutoff(profile)
        full_discharged_as = cumulative_trapezoid(full_current_a, full_time_s, initial=0.0)
        predicted_soc_full = 1.0 - full_discharged_as / (3600.0 * predicted_capacity_ah)
        predicted_voltage_full = relation(predicted_soc_full)
        full_error_mv = (predicted_voltage_full - measured_voltage_v) * 1000.0
        full_support_errors.append(full_error_mv)
        full_metrics = error_metrics(full_error_mv)

        # Primary common support ends at the earlier measured/predicted endpoint.
        common_end_s = profile.measured_cutoff_s
        if math.isfinite(predicted_endpoint_s):
            common_end_s = min(common_end_s, predicted_endpoint_s)
        common_time_s, common_current_a, common_measured_voltage_v = arrays_to_time(
            profile, common_end_s
        )
        common_discharged_as = cumulative_trapezoid(
            common_current_a, common_time_s, initial=0.0
        )
        predicted_soc_common = 1.0 - common_discharged_as / (
            3600.0 * predicted_capacity_ah
        )
        predicted_voltage_common = relation(predicted_soc_common)
        common_error_mv = (predicted_voltage_common - common_measured_voltage_v) * 1000.0
        common_support_errors.append(common_error_mv)
        common_metrics = error_metrics(common_error_mv)

        # Evaluation-only diagnostic: isolate U_dev mapping error by substituting
        # the observed cycle capacity. This quantity is never used for prediction.
        observed_capacity_soc = 1.0 - common_discharged_as / (
            3600.0 * profile.capacity_to_cutoff_ah
        )
        observed_capacity_voltage = relation(observed_capacity_soc)
        observed_capacity_error_mv = (
            observed_capacity_voltage - common_measured_voltage_v
        ) * 1000.0
        observed_capacity_diagnostic_errors.append(observed_capacity_error_mv)
        observed_capacity_metrics = error_metrics(observed_capacity_error_mv)

        endpoint_error_min = (
            (predicted_endpoint_s - profile.measured_cutoff_s) / 60.0
            if math.isfinite(predicted_endpoint_s)
            else math.nan
        )
        common_support_fraction = (
            common_end_s / profile.measured_cutoff_s
            if profile.measured_cutoff_s > 0.0
            else math.nan
        )

        cycle_rows.append(
            {
                "window": window_name,
                "cycle": profile.cycle,
                "predicted_capacity_mAh": predicted_capacity_ah * 1000.0,
                "observed_capacity_to_3V_mAh_evaluation_only": profile.capacity_to_cutoff_ah * 1000.0,
                "common_support_end_min": common_end_s / 60.0,
                "common_support_time_fraction": common_support_fraction,
                "voltage_MAE_mV": common_metrics["MAE"],
                "voltage_RMSE_mV": common_metrics["RMSE"],
                "voltage_bias_mV": common_metrics["bias"],
                "voltage_max_abs_mV": common_metrics["max_abs"],
                "full_measured_support_voltage_MAE_mV": full_metrics["MAE"],
                "full_measured_support_voltage_RMSE_mV": full_metrics["RMSE"],
                "full_measured_support_voltage_bias_mV": full_metrics["bias"],
                "full_measured_support_voltage_max_abs_mV": full_metrics["max_abs"],
                "diagnostic_observed_capacity_voltage_MAE_mV": observed_capacity_metrics["MAE"],
                "diagnostic_observed_capacity_voltage_RMSE_mV": observed_capacity_metrics["RMSE"],
                "diagnostic_observed_capacity_voltage_bias_mV": observed_capacity_metrics["bias"],
                "observed_endpoint_time_min": profile.measured_cutoff_s / 60.0,
                "predicted_capacity_derived_endpoint_time_min": (
                    predicted_endpoint_s / 60.0
                    if math.isfinite(predicted_endpoint_s)
                    else math.nan
                ),
                "endpoint_time_error_min": endpoint_error_min,
                "predicted_empirical_voltage_crossing_min": (
                    t_voltage / 60.0 if math.isfinite(t_voltage) else math.nan
                ),
                "predicted_SOC_depletion_min": (
                    t_soc / 60.0 if math.isfinite(t_soc) else math.nan
                ),
                "endpoint_event_type": event_type,
                "endpoint_evaluable": bool(math.isfinite(predicted_endpoint_s)),
                "trajectory_points_full_measured_support": int(len(full_time_s)),
                "trajectory_points_primary_common_support": int(len(common_time_s)),
            }
        )

        # Save the full measured support for transparent plotting. A flag marks
        # points included in the primary common-support metrics.
        observed_capacity_soc_full = 1.0 - full_discharged_as / (
            3600.0 * profile.capacity_to_cutoff_ah
        )
        observed_capacity_voltage_full = relation(observed_capacity_soc_full)
        observed_capacity_error_full_mv = (
            observed_capacity_voltage_full - measured_voltage_v
        ) * 1000.0
        for index in range(len(full_time_s)):
            soc_value = float(predicted_soc_full[index])
            soc_bin = pd.cut(
                pd.Series([soc_value]),
                bins=SOC_BINS,
                labels=SOC_LABELS,
                include_lowest=True,
                right=False,
            ).iloc[0]
            long_rows.append(
                {
                    "window": window_name,
                    "cycle": profile.cycle,
                    "time_s": float(full_time_s[index]),
                    "current_a": float(full_current_a[index]),
                    "measured_voltage_v": float(measured_voltage_v[index]),
                    "predicted_voltage_v": float(predicted_voltage_full[index]),
                    "predicted_soc": soc_value,
                    "soc_region": str(soc_bin) if not pd.isna(soc_bin) else "outside_0_1",
                    "voltage_error_mV": float(full_error_mv[index]),
                    "within_primary_common_support": bool(
                        full_time_s[index] <= common_end_s + 1e-9
                    ),
                    "evaluation_only_observed_capacity_soc": float(
                        observed_capacity_soc_full[index]
                    ),
                    "evaluation_only_observed_capacity_voltage_v": float(
                        observed_capacity_voltage_full[index]
                    ),
                    "evaluation_only_observed_capacity_error_mV": float(
                        observed_capacity_error_full_mv[index]
                    ),
                }
            )

    cycle_df = pd.DataFrame(cycle_rows)
    long_df = pd.DataFrame(long_rows)
    pooled_common = error_metrics(np.concatenate(common_support_errors))
    pooled_full = error_metrics(np.concatenate(full_support_errors))
    pooled_observed_capacity = error_metrics(
        np.concatenate(observed_capacity_diagnostic_errors)
    )

    endpoint_error = cycle_df["endpoint_time_error_min"].to_numpy(dtype=float)
    finite_endpoint = endpoint_error[np.isfinite(endpoint_error)]

    summary_df = pd.DataFrame(
        [
            {
                "window": window_name,
                "cycles": int(len(cycle_df)),
                "primary_voltage_support": "common support to earlier observed/predicted endpoint",
                "common_support_voltage_MAE_mV": pooled_common["MAE"],
                "common_support_voltage_RMSE_mV": pooled_common["RMSE"],
                "common_support_voltage_bias_mV": pooled_common["bias"],
                "common_support_voltage_max_abs_mV": pooled_common["max_abs"],
                "cycle_mean_common_support_voltage_MAE_mV": float(
                    cycle_df["voltage_MAE_mV"].mean()
                ),
                "cycle_mean_common_support_voltage_RMSE_mV": float(
                    cycle_df["voltage_RMSE_mV"].mean()
                ),
                "mean_common_support_time_fraction": float(
                    cycle_df["common_support_time_fraction"].mean()
                ),
                "full_measured_support_voltage_MAE_mV_diagnostic": pooled_full["MAE"],
                "full_measured_support_voltage_RMSE_mV_diagnostic": pooled_full["RMSE"],
                "full_measured_support_voltage_bias_mV_diagnostic": pooled_full["bias"],
                "full_measured_support_voltage_max_abs_mV_diagnostic": pooled_full["max_abs"],
                "observed_capacity_diagnostic_voltage_MAE_mV": pooled_observed_capacity["MAE"],
                "observed_capacity_diagnostic_voltage_RMSE_mV": pooled_observed_capacity["RMSE"],
                "observed_capacity_diagnostic_voltage_bias_mV": pooled_observed_capacity["bias"],
                "endpoint_evaluable_count": int(np.sum(cycle_df["endpoint_evaluable"])),
                "endpoint_time_MAE_min": (
                    float(np.mean(np.abs(finite_endpoint)))
                    if len(finite_endpoint)
                    else math.nan
                ),
                "endpoint_time_RMSE_min": (
                    float(np.sqrt(np.mean(finite_endpoint**2)))
                    if len(finite_endpoint)
                    else math.nan
                ),
                "endpoint_time_bias_min": (
                    float(np.mean(finite_endpoint)) if len(finite_endpoint) else math.nan
                ),
                "SOC_first_count": int(
                    np.sum(cycle_df["endpoint_event_type"] == "SOC-first")
                ),
                "voltage_first_count": int(
                    np.sum(cycle_df["endpoint_event_type"] == "voltage-first")
                ),
                "tie_count_within_1s": int(
                    np.sum(cycle_df["endpoint_event_type"] == "tie")
                ),
                "not_reached_count": int(
                    np.sum(cycle_df["endpoint_event_type"] == "not-reached")
                ),
                "voltage_bias_definition": "prediction minus measurement",
                "endpoint_interpretation": (
                    "capacity-derived endpoint under U_dev(0)=3.0 V construction; "
                    "not independent voltage-triggered cutoff validation"
                ),
            }
        ]
    )

    soc_rows: list[dict[str, Any]] = []
    primary_long = long_df[long_df["within_primary_common_support"]].copy()
    for region in SOC_LABELS:
        subset = primary_long[primary_long["soc_region"] == region]
        region_metrics = error_metrics(subset["voltage_error_mV"].to_numpy(dtype=float))
        soc_rows.append(
            {
                "window": window_name,
                "SOC_region": region,
                "support": "primary common support",
                "samples": int(len(subset)),
                "voltage_MAE_mV": region_metrics["MAE"],
                "voltage_RMSE_mV": region_metrics["RMSE"],
                "voltage_bias_mV": region_metrics["bias"],
                "voltage_max_abs_mV": region_metrics["max_abs"],
            }
        )
    soc_df = pd.DataFrame(soc_rows)
    return cycle_df, long_df, summary_df, soc_df

def moving_block_bootstrap_mean_ci(
    values: np.ndarray,
    *,
    samples: int,
    block_length: int,
    seed: int,
) -> dict[str, float | int | str]:
    """Circular moving-block bootstrap CI for an ordered cycle-level mean."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {
            "mean": math.nan,
            "ci95_low": math.nan,
            "ci95_high": math.nan,
            "n": 0,
            "method": "circular moving-block bootstrap",
            "block_length": 0,
        }
    n = len(x)
    block = max(1, min(int(block_length), n))
    blocks_per_sample = int(math.ceil(n / block))
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, n, size=(samples, blocks_per_sample))
    offsets = np.arange(block, dtype=int)
    indices = (starts[..., None] + offsets) % n
    indices = indices.reshape(samples, -1)[:, :n]
    means = x[indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return {
        "mean": float(np.mean(x)),
        "ci95_low": float(low),
        "ci95_high": float(high),
        "n": int(n),
        "method": "circular moving-block bootstrap",
        "block_length": int(block),
    }


def build_max_error_diagnostics(predictions: pd.DataFrame) -> pd.DataFrame:
    """Locate the largest absolute residual on primary and full supports."""
    rows: list[dict[str, Any]] = []
    for window in predictions["window"].drop_duplicates():
        window_df = predictions[predictions["window"] == window]
        for support_name, subset in [
            (
                "primary_common_support",
                window_df[window_df["within_primary_common_support"]],
            ),
            ("full_measured_support_diagnostic", window_df),
        ]:
            if subset.empty:
                continue
            index = subset["voltage_error_mV"].abs().idxmax()
            row = subset.loc[index]
            rows.append(
                {
                    "window": window,
                    "support": support_name,
                    "cycle": int(row["cycle"]),
                    "time_min": float(row["time_s"] / 60.0),
                    "predicted_SOC": float(row["predicted_soc"]),
                    "measured_voltage_V": float(row["measured_voltage_v"]),
                    "predicted_voltage_V": float(row["predicted_voltage_v"]),
                    "prediction_minus_measurement_error_mV": float(
                        row["voltage_error_mV"]
                    ),
                    "absolute_error_mV": float(abs(row["voltage_error_mV"])),
                    "at_discharge_onset": bool(abs(row["time_s"]) <= 1e-9),
                }
            )
    return pd.DataFrame(rows)

def configure_plotting() -> None:
    font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    if Path(font_path).is_file():
        font_manager.fontManager.addfont(font_path)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 9.5,
            "axes.labelsize": 10,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def draw_main_figure(
    final_long: pd.DataFrame,
    final_cycle: pd.DataFrame,
    final_soc: pd.DataFrame,
    representative_cycle: int,
    output_png: Path,
    output_pdf: Path,
    cell_label: str,
) -> None:
    configure_plotting()
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2))
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        bottom=0.10,
        top=0.91,
        wspace=0.30,
        hspace=0.38,
    )

    representative = final_long[final_long["cycle"] == representative_cycle].copy()
    representative_metric = final_cycle[
        final_cycle["cycle"] == representative_cycle
    ].iloc[0]
    measured_endpoint = float(representative_metric["observed_endpoint_time_min"])
    predicted_endpoint = float(
        representative_metric["predicted_capacity_derived_endpoint_time_min"]
    )

    ax = axes[0, 0]
    ax.plot(
        representative["time_s"] / 60.0,
        representative["measured_voltage_v"],
        linewidth=1.8,
        label="Measured voltage",
    )
    ax.plot(
        representative["time_s"] / 60.0,
        representative["predicted_voltage_v"],
        linewidth=1.5,
        linestyle="--",
        label="Empirical reconstruction",
    )
    ax.axhline(V_CUT, linestyle=":", linewidth=1.0, label="3.0-V boundary")
    ax.set_xlabel("Discharge time (min)")
    ax.set_ylabel("Terminal voltage (V)")
    ax.set_title(
        f"Representative final-confirmation cycle: Cycle {representative_cycle}",
        
    )
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(
        representative["time_s"] / 60.0,
        representative["voltage_error_mV"],
        linewidth=1.3,
        color="tab:blue",
    )
    ax.axhline(0.0, color="0.45", linestyle="--", linewidth=0.9)
    if math.isfinite(predicted_endpoint):
        ax.axvline(
            predicted_endpoint,
            color="tab:orange",
            linestyle="--",
            linewidth=1.4,
            label="Predicted capacity-derived endpoint",
        )
    ax.axvline(
        measured_endpoint,
        color="tab:blue",
        linestyle=":",
        linewidth=1.6,
        label="Measured 3.0-V endpoint",
    )
    if math.isfinite(predicted_endpoint) and predicted_endpoint < measured_endpoint:
        ax.axvspan(
            predicted_endpoint,
            measured_endpoint,
            color="0.85",
            alpha=0.35,
            label="Outside primary common support",
        )
    ax.set_xlabel("Discharge time (min)")
    ax.set_ylabel("Prediction minus measurement residual (mV)")
    ax.set_title("(b) Voltage residual and endpoint locations")
    ax.legend(frameon=False, fontsize=7.4)

    ax = axes[1, 0]
    x = np.arange(len(final_soc))
    width = 0.36
    ax.bar(x - width / 2, final_soc["voltage_MAE_mV"], width, label="MAE")
    ax.bar(x + width / 2, final_soc["voltage_RMSE_mV"], width, label="RMSE")
    ax.set_xticks(x)
    ax.set_xticklabels(final_soc["SOC_region"])
    ax.set_xlabel("Predicted SOC region")
    ax.set_ylabel("Voltage error on common support (mV)")
    ax.set_title("(c) Voltage error across SOC regions")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    endpoint_df = final_cycle[
        final_cycle["predicted_capacity_derived_endpoint_time_min"].notna()
    ].copy()
    if endpoint_df.empty:
        ax.text(
            0.5,
            0.5,
            "No evaluable predicted endpoint in this window",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
    else:
        observed_time = endpoint_df["observed_endpoint_time_min"].to_numpy(float)
        predicted_time = endpoint_df[
            "predicted_capacity_derived_endpoint_time_min"
        ].to_numpy(float)
        endpoint_error = predicted_time - observed_time
        endpoint_mae = float(np.mean(np.abs(endpoint_error)))
        endpoint_bias = float(np.mean(endpoint_error))

        ax.scatter(observed_time, predicted_time, s=30, color="tab:blue")
        lower = float(min(observed_time.min(), predicted_time.min()))
        upper = float(max(observed_time.max(), predicted_time.max()))
        pad = max(0.10, 0.05 * (upper - lower))
        lower -= pad
        upper += pad
        limits = [lower, upper]
        ax.plot(limits, limits, color="0.40", linestyle="--", linewidth=1.0)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_aspect("equal", adjustable="box")
        ax.text(
            0.04,
            0.96,
            f"MAE = {endpoint_mae:.3f} min\nBias = {endpoint_bias:.3f} min",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8.0,
        )
    ax.set_xlabel("Measured 3.0-V endpoint time (min)")
    ax.set_ylabel("Predicted capacity-derived endpoint time (min)")
    ax.set_title("(d) Capacity-derived endpoint-time agreement")

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linewidth=0.5, alpha=0.35)

    fig.suptitle(
        f"Confirmatory Evaluation of Empirical Terminal-Voltage Reconstruction for {cell_label}",
        fontsize=13,
    )
    fig.savefig(output_png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def draw_relation_figure(
    relation_fit: EmpiricalRelation,
    relation_outer: EmpiricalRelation,
    output_png: Path,
    output_pdf: Path,
) -> None:
    configure_plotting()
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.plot(
        relation_fit.soc_grid,
        relation_fit.voltage_grid,
        linestyle="--",
        linewidth=1.5,
        label="Cycles 1-40 relation",
    )
    ax.plot(
        relation_outer.soc_grid,
        relation_outer.voltage_grid,
        linewidth=1.8,
        label="Cycles 1-60 frozen relation",
    )
    ax.axhline(V_CUT, linestyle=":", linewidth=1.0, label="3.0-V construction boundary")
    ax.set_xlabel("SOC")
    ax.set_ylabel("Terminal voltage (V)")
    ax.set_title("Development-Derived Loaded Terminal-Voltage-SOC Relations")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cell-label", type=str, default="CS2-21")
    args = parser.parse_args()

    raw_file = args.raw_file.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if not raw_file.is_file():
        raise FileNotFoundError(raw_file)

    profiles, extraction_audit = read_profiles(raw_file)
    cycles = sorted(profiles)
    fit_1_40 = [profiles[c] for c in cycles if c <= FIT_END]
    validation_41_60 = [profiles[c] for c in cycles if FIT_END < c <= VALIDATION_END]
    outer_1_60 = [profiles[c] for c in cycles if c <= VALIDATION_END]
    robustness_61_80 = [profiles[c] for c in cycles if VALIDATION_END < c <= ROBUSTNESS_END]
    final_81_101 = [profiles[c] for c in cycles if c > ROBUSTNESS_END]
    if not all([fit_1_40, validation_41_60, outer_1_60, robustness_61_80, final_81_101]):
        raise ExperimentError("One or more chronological windows are empty.")

    capacity_fit = fit_saturating_capacity(fit_1_40)
    relation_fit = build_empirical_relation(fit_1_40)
    capacity_outer = fit_saturating_capacity(outer_1_60)
    relation_outer = build_empirical_relation(outer_1_60)

    window_results = []
    for name, evaluation_profiles, capacity, relation in [
        ("validation_41_60", validation_41_60, capacity_fit, relation_fit),
        ("robustness_61_80", robustness_61_80, capacity_outer, relation_outer),
        ("final_confirmation_81_101", final_81_101, capacity_outer, relation_outer),
    ]:
        window_results.append(
            evaluate_window(name, evaluation_profiles, capacity, relation)
        )

    cycle_metrics = pd.concat([result[0] for result in window_results], ignore_index=True)
    predictions = pd.concat([result[1] for result in window_results], ignore_index=True)
    summary = pd.concat([result[2] for result in window_results], ignore_index=True)
    soc_metrics = pd.concat([result[3] for result in window_results], ignore_index=True)

    relation_frame = pd.DataFrame(
        {
            "soc": relation_outer.soc_grid,
            "U_dev_voltage_v": relation_outer.voltage_grid,
        }
    )
    parameter_frame = pd.DataFrame(
        [
            {"fit_scope": "cycles_1_40", **capacity_fit, "U_dev_at_SOC0_V": float(relation_fit(0.0)), "U_dev_at_SOC1_V": float(relation_fit(1.0)), "relation_source": relation_fit.source},
            {"fit_scope": "cycles_1_60_frozen", **capacity_outer, "U_dev_at_SOC0_V": float(relation_outer(0.0)), "U_dev_at_SOC1_V": float(relation_outer(1.0)), "relation_source": relation_outer.source},
        ]
    )

    final_cycle = cycle_metrics[cycle_metrics["window"] == "final_confirmation_81_101"].copy()
    final_long = predictions[predictions["window"] == "final_confirmation_81_101"].copy()
    final_soc = soc_metrics[soc_metrics["window"] == "final_confirmation_81_101"].copy()

    median_mae = float(final_cycle["voltage_MAE_mV"].median())
    representative_cycle = int(
        final_cycle.iloc[(final_cycle["voltage_MAE_mV"] - median_mae).abs().argmin()]["cycle"]
    )

    confidence_rows = []
    for metric_name, values in [
        (
            "cycle_mean_common_support_voltage_MAE_mV",
            final_cycle["voltage_MAE_mV"].to_numpy(float),
        ),
        (
            "cycle_mean_common_support_voltage_RMSE_mV",
            final_cycle["voltage_RMSE_mV"].to_numpy(float),
        ),
        (
            "capacity_derived_endpoint_absolute_error_min",
            np.abs(final_cycle["endpoint_time_error_min"].to_numpy(float)),
        ),
    ]:
        ci = moving_block_bootstrap_mean_ci(
            values,
            samples=BOOTSTRAP_SAMPLES,
            block_length=MOVING_BLOCK_LENGTH,
            seed=BOOTSTRAP_SEED,
        )
        confidence_rows.append(
            {
                "window": "final_confirmation_81_101",
                "metric": metric_name,
                **ci,
                "bootstrap_samples": BOOTSTRAP_SAMPLES,
                "seed": BOOTSTRAP_SEED,
            }
        )
    confidence = pd.DataFrame(confidence_rows)

    max_error_diagnostics = build_max_error_diagnostics(predictions)
    final_summary = summary[
        summary["window"] == "final_confirmation_81_101"
    ].iloc[0]
    capacity_mapping_diagnostic = pd.DataFrame(
        [
            {
                "setting": "Frozen capacity prediction + U_dev",
                "role": "primary end-to-end reconstruction",
                "voltage_MAE_mV": final_summary[
                    "common_support_voltage_MAE_mV"
                ],
                "voltage_RMSE_mV": final_summary[
                    "common_support_voltage_RMSE_mV"
                ],
                "voltage_bias_mV": final_summary[
                    "common_support_voltage_bias_mV"
                ],
            },
            {
                "setting": "Observed capacity + U_dev",
                "role": "evaluation-only diagnostic; observed capacity is not a model input",
                "voltage_MAE_mV": final_summary[
                    "observed_capacity_diagnostic_voltage_MAE_mV"
                ],
                "voltage_RMSE_mV": final_summary[
                    "observed_capacity_diagnostic_voltage_RMSE_mV"
                ],
                "voltage_bias_mV": final_summary[
                    "observed_capacity_diagnostic_voltage_bias_mV"
                ],
            },
        ]
    )

    summary.to_csv(output_dir / "Table3_terminal_voltage_reconstruction_summary.csv", index=False)
    # Legacy filename retained for manuscript workflows, but content now follows
    # the corrected endpoint-consistency interpretation.
    summary.to_csv(output_dir / "Table3_empirical_voltage_summary.csv", index=False)
    soc_metrics.to_csv(output_dir / "TableS_terminal_SOC_error.csv", index=False)
    cycle_metrics.to_csv(output_dir / "empirical_voltage_metrics_by_cycle.csv", index=False)
    predictions.to_csv(output_dir / "empirical_voltage_predictions_long.csv", index=False)
    relation_frame.to_csv(output_dir / "frozen_Udev_SOC_relation.csv", index=False)
    parameter_frame.to_csv(output_dir / "empirical_voltage_model_parameters.csv", index=False)
    confidence.to_csv(output_dir / "final_confirmation_moving_block_bootstrap_CI.csv", index=False)
    max_error_diagnostics.to_csv(output_dir / "max_error_diagnostics.csv", index=False)
    capacity_mapping_diagnostic.to_csv(
        output_dir / "TableS_capacity_mapping_diagnostic.csv", index=False
    )

    draw_main_figure(
        final_long,
        final_cycle,
        final_soc,
        representative_cycle,
        output_dir / "Fig_empirical_voltage_final_confirmation.png",
        output_dir / "Fig_empirical_voltage_final_confirmation.pdf",
        args.cell_label,
    )
    draw_relation_figure(
        relation_fit,
        relation_outer,
        output_dir / "Fig_Udev_SOC_relation.png",
        output_dir / "Fig_Udev_SOC_relation.pdf",
    )

    audit = {
        "source": {
            "cell_label": args.cell_label,
            "raw_file": str(raw_file),
            "raw_file_bytes": raw_file.stat().st_size,
            "raw_file_sha256": sha256_file(raw_file),
        },
        "extraction": extraction_audit,
        "chronological_protocol": {
            "fit_cycles": [profile.cycle for profile in fit_1_40],
            "validation_cycles": [profile.cycle for profile in validation_41_60],
            "frozen_refit_cycles": [profile.cycle for profile in outer_1_60],
            "robustness_cycles": [profile.cycle for profile in robustness_61_80],
            "final_confirmation_cycles": [profile.cycle for profile in final_81_101],
            "test_measured_capacity_used_as_model_input": False,
            "test_measured_voltage_used_for_fitting": False,
        },
        "model": {
            "formula": "V_hat(t) = U_dev(z_hat(t))",
            "SOC_formula": "z_hat(t)=1-integral(I dt)/(3600*Q_hat(N))",
            "U_dev_definition": relation_outer.source,
            "not_claimed_as": ["OCV", "equilibrium voltage", "independently calibrated OCV curve"],
            "electrical_parameters_removed": ["R0", "R1", "C1", "v1"],
            "voltage_layer_role": "held-out terminal-voltage reconstruction and endpoint-consistency audit, not a validated dynamic equivalent-circuit model",
            "primary_voltage_error_support": "common support ending at the earlier measured/predicted endpoint",
            "voltage_bias_definition": "prediction minus measurement",
            "endpoint_interpretation": "capacity-derived endpoint because U_dev(0)=3.0 V is imposed by construction",
            "V_cut_V": V_CUT,
            "SOC_min": SOC_MIN,
            "tie_tolerance_s": TIE_TOLERANCE_S,
        },
        "parameters": parameter_frame.to_dict(orient="records"),
        "summary_metrics": summary.to_dict(orient="records"),
        "terminal_SOC_metrics": soc_metrics.to_dict(orient="records"),
        "moving_block_bootstrap_confidence_intervals": confidence.to_dict(orient="records"),
        "capacity_mapping_diagnostic": capacity_mapping_diagnostic.to_dict(orient="records"),
        "maximum_error_diagnostics": max_error_diagnostics.to_dict(orient="records"),
        "representative_cycle_selection": {
            "rule": "final-confirmation cycle whose per-cycle voltage MAE is closest to the final-window median",
            "selected_cycle": representative_cycle,
        },
        "final_claim": (
            "The evidence supports chronological capacity aging plus held-out reconstruction of "
            "loaded terminal-voltage trajectories. Because U_dev(0)=3.0 V is imposed by construction, "
            "the reported endpoint timing is capacity-derived and voltage/SOC crossings are mathematically "
            "coupled. The voltage layer is therefore an endpoint-consistency audit, not independent "
            "validation of voltage-triggered cutoff or load-induced voltage sag."
        ),
    }
    with (output_dir / "empirical_voltage_audit.json").open("w", encoding="utf-8") as stream:
        json.dump(json_ready(audit), stream, ensure_ascii=False, indent=2)

    print("Route-A terminal-voltage reconstruction and endpoint-consistency audit completed.")
    print(f"Raw SHA-256: {audit['source']['raw_file_sha256']}")
    print(f"Accepted cycles: {extraction_audit['accepted_cycles']}")
    print(f"Final confirmation cycles: {[profile.cycle for profile in final_81_101]}")
    print(summary.to_string(index=False))
    print(f"Representative cycle: {representative_cycle}")
    print(f"Results written to: {output_dir}")


if __name__ == "__main__":
    main()
