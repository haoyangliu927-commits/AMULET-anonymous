from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd


@dataclass
class GeneAmbiguityResult:
    table: pd.DataFrame
    ambiguity_weight: np.ndarray
    raw_ambiguity_score: np.ndarray
    robust_z: np.ndarray


def _as_bool_matrix(matrix, *, name: str, shape: tuple[int, int] | None = None) -> np.ndarray:
    values = np.asarray(matrix, dtype=bool)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError(f"{name} must be a square matrix, got shape {values.shape}.")
    if shape is not None and values.shape != shape:
        raise ValueError(f"{name} shape {values.shape} does not match expected {shape}.")
    return values


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60.0, 60.0)))


def _rank_normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.size <= 1:
        return np.zeros_like(values, dtype=np.float64)
    return pd.Series(values).rank(method="average", pct=True).to_numpy(dtype=np.float64)


def compute_gene_ambiguity(
    positive_mask,
    negative_mask,
    partial_pos_mask=None,
    directed_inclusion_mask=None,
    gene_names: Sequence[str] | None = None,
    *,
    robust_center_z: float = 3.0,
    robust_temperature: float = 0.75,
    contradiction_weight: float = 1.0,
    bidirectional_inclusion_weight: float = 0.3,
    eps: float = 1e-8,
) -> GeneAmbiguityResult:
    """Score genes with contradictory rather than merely diverse relation evidence."""
    positive = _as_bool_matrix(positive_mask, name="positive_mask")
    negative = _as_bool_matrix(negative_mask, name="negative_mask", shape=positive.shape)
    partial = (
        np.zeros_like(positive, dtype=bool)
        if partial_pos_mask is None
        else _as_bool_matrix(partial_pos_mask, name="partial_pos_mask", shape=positive.shape)
    )
    directed = (
        np.zeros_like(positive, dtype=bool)
        if directed_inclusion_mask is None
        else _as_bool_matrix(
            directed_inclusion_mask,
            name="directed_inclusion_mask",
            shape=positive.shape,
        )
    )

    np.fill_diagonal(positive, False)
    np.fill_diagonal(negative, False)
    np.fill_diagonal(partial, False)
    np.fill_diagonal(directed, False)

    contradiction_mask_pos_neg = positive & negative
    contradiction_mask_partial_neg = partial & negative
    contradiction_degree = (
        contradiction_mask_pos_neg.sum(axis=1) + contradiction_mask_partial_neg.sum(axis=1)
    ).astype(np.float64)
    contradiction_score = np.log1p(contradiction_degree)

    bidirectional_inclusion = directed & directed.T
    bidirectional_inclusion_degree = bidirectional_inclusion.sum(axis=1).astype(np.float64)
    bidirectional_inclusion_score = np.log1p(bidirectional_inclusion_degree)

    contradiction_score_norm = _rank_normalize(contradiction_score)
    bidirectional_score_norm = _rank_normalize(bidirectional_inclusion_score)
    raw_score = (
        float(contradiction_weight) * contradiction_score_norm
        + float(bidirectional_inclusion_weight) * bidirectional_score_norm
    )

    median = float(np.median(raw_score))
    mad = float(np.median(np.abs(raw_score - median)))
    robust_z = (raw_score - median) / (1.4826 * mad + eps)
    ambiguity_weight = _sigmoid((robust_z - robust_center_z) / robust_temperature)

    if gene_names is None:
        names = [str(idx) for idx in range(positive.shape[0])]
    else:
        names = [str(name) for name in list(gene_names)[: positive.shape[0]]]

    table = pd.DataFrame(
        {
            "gene_index": np.arange(positive.shape[0], dtype=int),
            "gene_name": names,
            "contradiction_degree": contradiction_degree.astype(int),
            "contradiction_score": contradiction_score,
            "contradiction_score_norm": contradiction_score_norm,
            "bidirectional_inclusion_degree": bidirectional_inclusion_degree.astype(int),
            "bidirectional_inclusion_score": bidirectional_inclusion_score,
            "bidirectional_inclusion_score_norm": bidirectional_score_norm,
            "contradiction_weight": float(contradiction_weight),
            "bidirectional_inclusion_weight": float(bidirectional_inclusion_weight),
            "raw_ambiguity_score": raw_score,
            "robust_z": robust_z,
            "ambiguity_weight": ambiguity_weight,
        }
    )
    return GeneAmbiguityResult(
        table=table,
        ambiguity_weight=ambiguity_weight.astype(np.float32),
        raw_ambiguity_score=raw_score.astype(np.float32),
        robust_z=robust_z.astype(np.float32),
    )
