"""
rpy2_bridge.py — shared utilities for pandas <-> R matrix conversion.
Import this module in all Row scripts; it activates converters once at import time.
"""
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
import pandas as pd

pandas2ri.activate()
numpy2ri.activate()

base = importr("base")
immunedeconv = importr("immunedeconv")


def df_to_r_matrix(df: pd.DataFrame):
    """Convert genes×samples DataFrame to an R numeric matrix with rownames/colnames."""
    r_df = pandas2ri.py2rpy(df)
    r_mat = base.as_matrix(r_df)
    base.rownames(r_mat)
    base.colnames(r_mat)
    return r_mat


def r_to_df(r_obj) -> pd.DataFrame:
    """Convert an R matrix or data.frame result back to a pandas DataFrame."""
    return pandas2ri.rpy2py(r_obj)
