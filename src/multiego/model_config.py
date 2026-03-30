# model_config.py
from dataclasses import dataclass


@dataclass(frozen=True)
class ModelConfig:
    """
    Configuration parameters for the interaction model.

    Attributes:
        max_bond_separation: Maximum number of bonds for applying the nt-bonds-rule beyond the typical 1-4 (3 bonds)
            (e.g., exclude interactions for atoms separated by <= this number).
    """

    max_bond_separation: int = 5
    bond14_separation: int = 3

    def __post_init__(self):
        assert self.max_bond_separation >= self.bond14_separation, (
            f"max_bond_separation ({self.max_bond_separation}) must be >= "
            f"bond14_separation ({self.bond14_separation})"
        )


config = ModelConfig()
