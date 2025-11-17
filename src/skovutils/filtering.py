"""
Filtering utilities for scikit-bio tables.
"""

from __future__ import annotations


# --- Minimum abundance and prevalence filtering ---
def minmin(table, abundance: float = 0.01, prevalence: float = 0.1):
    """
    Filter a skbio Table to only include observations that meet minimum
    abundance and prevalence thresholds.

    Parameters
    ----------
    table : skbio.Table
        Table to be filtered (observations x samples).
    abundance : float, default=0.01
        Minimum abundance threshold for inclusion.
    prevalence : float, default=0.1
        Fraction of samples that must exceed ``abundance``.

    Returns
    -------
    skbio.Table
        Filtered table containing only the observations that meet the
        abundance and prevalence criteria.
    """
    # Convert table to DataFrame
    df = table.to_dataframe()

    # Create boolean mask where abundance is above threshold
    abundance_mask = df >= abundance

    # Convert boolean mask to integers (True = 1, False = 0)
    abundance_mask = abundance_mask.astype(int)

    # Count how many samples meet the abundance threshold per observation
    prevalence_counts = abundance_mask.sum(axis=1)

    # Calculate the minimum number of samples required for prevalence
    prevalence_threshold = prevalence * df.shape[1]

    # Identify observations that meet the prevalence threshold
    observations_to_keep = prevalence_counts[
        prevalence_counts >= prevalence_threshold
    ].index.tolist()

    # Filter the table to include only selected observations
    table_filtered = table.filter(observations_to_keep, axis="observation")

    return table_filtered


if __name__ == "__main__":  # pragma: no cover - usage example
    from skbio import Table
    import numpy as np

    # Create a sample skbio Table
    data = np.array(
        [
            [0.0, 0.02, 0.03],
            [0.1, 0.0, 0.0],
            [0.05, 0.06, 0.07],
            [0.0, 0.0, 0.0],
        ]
    )
    obs_ids = ["obs1", "obs2", "obs3", "obs4"]
    samp_ids = ["samp1", "samp2", "samp3"]
    table = Table(data, obs_ids, samp_ids)

    # Apply minmin filtering
    filtered_table = minmin(table, abundance=0.05, prevalence=0.5)

    # Print the filtered table
    print(filtered_table)
