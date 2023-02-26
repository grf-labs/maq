import numpy as np

from maq.ext import solver_cpp

class MAQ:
  r"""
  abc etc
  """

  def __init__(self, num_bootstrap=200, seed=42):
    self.num_bootstrap = num_bootstrap
    self.seed = seed

  def fit(reward, cost, budget):
    #etc

    return self
