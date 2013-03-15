import numpy as np
from symodesys.helpers import PieceWiseShiftedPolyTraj

def test_PieceWiseShiftedPolyTraj():
    pt1 = PieceWiseShiftedPolyTraj([0, 1.0], [[0, 0], [1, 3]])
    assert np.allclose(pt1._coeff_shifted_scaled, [0, 0, 0, 1])

def manual_PieceWiseShiftedPolyTraj():
    """
    Inspect manually
    """
    import matplotlib.pyplot as plt

    x = np.linspace(0, 1.0, 100)
    y = x ** 3
    pt1 = PieceWiseShiftedPolyTraj([0, 1.0], [[0, 0], [1, 3]])
    plt.plot(x, y)
    plt.plot(x, pt1(x))
