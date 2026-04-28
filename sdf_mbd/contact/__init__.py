"""SDF contact geometry, penalty baseline, and NCP utilities."""

from .ncp import complementarity_error, smooth_fischer_burmeister, smooth_fischer_burmeister_grad
from .penalty import PenaltyContactResult, penalty_normal_force_from_sdf
from .sdf_contact import finite_difference_gap_gradient, point_contact_jacobian, sdf_gap, sdf_normal

__all__ = [
    "PenaltyContactResult",
    "complementarity_error",
    "finite_difference_gap_gradient",
    "penalty_normal_force_from_sdf",
    "point_contact_jacobian",
    "sdf_gap",
    "sdf_normal",
    "smooth_fischer_burmeister",
    "smooth_fischer_burmeister_grad",
]

