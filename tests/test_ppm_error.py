import pandas as pd
from pathlib import Path


def test_ppm_error_matches_cli_output():
    data_file = Path(__file__).parent / "data" / "sage_cli_output.tsv"
    cli_df = pd.read_csv(data_file, sep="\t")
    calculated_ppm = ((cli_df["expmass"] - cli_df["calcmass"]) / cli_df["calcmass"]) * 1_000_000
    calculated_ppm.name = "ppm_error"
    pd.testing.assert_series_equal(calculated_ppm, cli_df["ppm_error"])
