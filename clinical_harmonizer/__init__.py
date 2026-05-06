from .base import ClinicalHarmonizer, SCHEMA_COLUMNS
from .tcga import TCGAHarmonizer

HARMONIZERS: dict[str, type[ClinicalHarmonizer]] = {
    "tcga": TCGAHarmonizer,
}

__all__ = ["ClinicalHarmonizer", "TCGAHarmonizer", "HARMONIZERS", "SCHEMA_COLUMNS"]
