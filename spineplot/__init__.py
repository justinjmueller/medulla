"""
Spineplot — plotting utilities for the Medulla analysis framework.
"""

__version__ = "1.0.4"
__author__ = "Justin Mueller"

# Expose key top-level classes
from .analysis import Analysis
from .core import (
    SpineFigure,
    SimpleFigure,
    Sample,
    Style,
    Variable,
)
from .artists import (
    SpineSpectra1D,
    SpineSpectra2D,
    ROCCurve,
    ConfusionMatrix,
    Ternary,
    SpineEfficiency,
)

__all__ = [
    "Analysis",
    "SpineFigure",
    "SimpleFigure",
    "SpineSpectra1D",
    "SpineSpectra2D",
    "ROCCurve",
    "ConfusionMatrix",
    "Ternary",
    "SpineEfficiency",
    "Style",
    "Variable",
    "Sample",
]

__all__ = [
    "Analysis",
    "SpineFigure",
    "SimpleFigure",
    "SpineSpectra1D",
    "SpineSpectra2D",
    "ROCCurve",
    "ConfusionMatrix",
    "Ternary",
    "SpineEfficiency",
    "Style",
    "Variable",
    "Sample",
]