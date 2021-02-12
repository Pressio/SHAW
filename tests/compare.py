
import numpy as np
import sys

if __name__== "__main__":
  print('Argument List:', str(sys.argv))

  file1 = str(sys.argv[1])
  file2 = str(sys.argv[2])
  tol   = float(sys.argv[3])
  firstLineHasSize = int(sys.argv[4])

  if firstLineHasSize:
    size1 = np.loadtxt(file1, max_rows=1)
    size2 = np.loadtxt(file2, max_rows=1)

    d1 = np.loadtxt(file1, skiprows=1)
    d2 = np.loadtxt(file2, skiprows=1)

    assert(np.isnan(d1).all() == False)
    assert(np.isnan(d2).all() == False)
    assert(np.allclose(d1, d2, atol=tol))
    assert(np.allclose(size1, size2))
    print("All good 1")

  else:
    d1 = np.loadtxt(file1, skiprows=1)
    d2 = np.loadtxt(file2, skiprows=1)
    assert(np.isnan(d1).all() == False)
    assert(np.isnan(d2).all() == False)
    assert(np.allclose(d1, d2, atol=tol))
    print("All good 2")
