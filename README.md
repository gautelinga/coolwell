# coolwell

This is the Git repository for the finite element implementation of the "cool well" project.

How to run:
```
mpiexec -n 4 python steady.py -fx 10.0 -kappa 0.01 -res 32
```

Dependencies:
* Python3
* Numpy
* FEniCS/Dolfin
* MeshPy
* h5py
