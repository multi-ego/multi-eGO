# model_config.py
from dataclasses import dataclass


@dataclass(frozen=True)
class ModelConfig:
    """
    Configuration parameters for the interaction model.

    Attributes:
        max_bond_separation: Maximum number of bonds for applying the nt-bonds-rule beyond the typical 1-4 (3 bonds)
            (e.g., exclude interactions for atoms separated by <= this number).
        bond14_separation: Number of bonds that defines a 1-4 interaction.
        epsilon_min: Minimum meaningful epsilon value (kJ/mol). Must be > 0.
            Used as the lower bound in the adaptive probability threshold formula
            and checked against per-reference epsilon values at startup.
        p_to_learn: Fraction of the total contact-probability mass that must be
            covered before the adaptive MD threshold is set. Should be close to 1
            (suggested: 0.9995). A warning is raised if it falls below 0.9.
        learn_tolerance: Relative deviation (unitless) between a trained LJ
            parameter and its MG prior below which the trained value is considered
            indistinguishable from the prior and is therefore discarded.
    """

    max_bond_separation: int = 5
    bond14_separation: int = 3
    epsilon_min: float = 0.07
    p_to_learn: float = 0.9995
    learn_tolerance: float = 0.01

    def __post_init__(self):
        assert self.max_bond_separation >= self.bond14_separation, (
            f"max_bond_separation ({self.max_bond_separation}) must be >= "
            f"bond14_separation ({self.bond14_separation})"
        )
        assert self.epsilon_min > 0.0, f"epsilon_min ({self.epsilon_min}) must be > 0"
        assert 0.0 < self.p_to_learn <= 1.0, f"p_to_learn ({self.p_to_learn}) must be in (0, 1]"
        if self.p_to_learn < 0.9:
            import warnings

            warnings.warn(
                f"p_to_learn ({self.p_to_learn}) is very small; suggested value is 0.9995",
                UserWarning,
                stacklevel=2,
            )


config = ModelConfig()
