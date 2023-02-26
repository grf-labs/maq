import pytest
import numpy as np
import numpy.testing as nt

from maq import MAQ
from maq.ext import solver_cpp

def test_tmp():
  r = np.random.rand(10, 3)
  c = np.random.rand(10, 3)
  solver_cpp(
    np.asfortranarray(r),
    np.asfortranarray(c),
    10,
    4,
    2,
    1
  )
