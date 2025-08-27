from collections import defaultdict
from typing import Dict, List, Union

from sagepy.core.scoring import Psm
from sagepy.ml import linear_discriminant, qvalue


def process_mzml(psm_collection: Union[List[Psm], Dict[str, List[Psm]]], prec_tol_ppm: float) -> List[Psm]:
    """Process a collection of PSMs and return those passing a 1% spectrum-level FDR.

    This function performs linear discriminant scoring and q-value estimation on a
    set of PSMs, keeps only the top-scoring PSM per spectrum group, and filters the
    results to a 1% spectrum-level false discovery rate.

    Args:
        psm_collection: Either a flat list of PSMs or a mapping of spectrum IDs to
            lists of PSM candidates.
        prec_tol_ppm: Precursor mass tolerance in parts-per-million, forwarded to the
            linear discriminant scoring routine.

    Returns:
        A list of PSMs passing the 1% spectrum-level FDR threshold.
    """

    # Collect all PSM objects into a single list
    all_features: List[Psm] = []
    if isinstance(psm_collection, dict):
        for psms in psm_collection.values():
            all_features.extend(psms)
    else:
        all_features.extend(psm_collection)

    # Score using the linear discriminant if available
    try:
        linear_discriminant.score_psms(all_features, prec_tol_ppm)
    except Exception:
        # If scoring fails (e.g., missing model), continue without modifying scores
        pass

    # Compute discriminant score and q-values for all PSMs
    qvalue.spectrum_q_value(all_features)

    # Perform spectrum-level FDR control
    grouped: Dict[str, List[Psm]] = defaultdict(list)
    for psm in all_features:
        grouped[getattr(psm, "_spec_group", psm.spec_idx)].append(psm)

    top_hits: List[Psm] = []
    for group_psms in grouped.values():
        group_psms.sort(key=lambda p: p.discriminant_score, reverse=True)
        top_hits.append(group_psms[0])

    # Filter by 1% spectrum-level FDR
    filtered = [psm for psm in top_hits if getattr(psm, "spectrum_q", 1.0) <= 0.01]

    return filtered
