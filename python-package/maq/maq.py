import numpy as np

from maq.ext import solver_cpp


class MAQ:
    r"""
    abc etc
    """

    def __init__(self, n_bootstrap=200, seed=42):
        self.n_bootstrap = n_bootstrap
        self.seed = seed

    def fit(reward, cost, budget):
        # etc

        return self
