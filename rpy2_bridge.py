"""
rpy2_bridge.py — shared utilities for pandas <-> R matrix conversion.

rpy2 >= 3.5 removed pandas2ri.activate() / numpy2ri.activate().
Conversion is now done explicitly inside localconverter() context blocks.
"""
import rpy2.robjects as ro
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import rpy2.robjects.pandas2ri as pandas2ri
import rpy2.robjects.numpy2ri as numpy2ri
import pandas as pd

# Combined converter — used inside localconverter() blocks below
converter = default_converter + pandas2ri.converter + numpy2ri.converter

base = importr("base")
immunedeconv = importr("immunedeconv")


def df_to_r_matrix(df: pd.DataFrame):
    """Convert a genes×samples DataFrame to an R numeric matrix."""
    with localconverter(converter):
        r_df = pandas2ri.py2rpy(df)
    r_mat = base.as_matrix(r_df)
    return r_mat


def r_to_df(r_obj) -> pd.DataFrame:
    """Convert an R matrix or data.frame result back to a pandas DataFrame."""
    with localconverter(converter):
        return pandas2ri.rpy2py(r_obj)
