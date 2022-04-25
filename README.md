# Coolwell

This is the Git repository for the finite element implementation of the "cool well" project.

Dependencies:
* Python3
* Numpy
* FEniCS/Dolfin
* MeshPy
* h5py

## Parameters

| parameter | description               | command line argument | default value | name in article[^1] |
|-----------|---------------------------|-----------------------|---------------|---------------------|
| res       | resolution <br>#cells ~ O((b+H)res/b) | `-res`                | 32            | res                 |
| L<sub>x     | width of the system including inlet and outlet | `-Lx`    | 2.0 | n/a |
| b         | width of the cavity <br>(inlet and outlet excluded)| `-b`   | 1.0 | h |
| H         | height of the channel | `-H` | 0.5 | w |
| f         |(horizontal) forcing      | `-fx`                 | 1.0           | f       |
|T<sub>top  | temperature of top surface | `-Ttop`    | 1.0 | T<sub>+ |
|T<sub>btm  | temperature of the bottom surface | `-Tbtm` | 0.0 | T<sub>- |
| &alpha;    | coefficient of thermal expansion | `-alpha` | 0.1 | &beta; |
| &nu;<sub>top | kinematic viscosity at temperature T<sub>top | `-nutop` | 1.0 | &nu; |
| &nu;<sub>btm | kinematic viscosity at temperature T<sub>btm | `-nubtm` | 2.0 | &nu; |
| k           | thermal diffusion coefficient | `-kappa` | 1.0 | k |
[^1]: This is the name used in [Baldelli et al. (2022)].

## How to run
```
mpiexec -np  [number of parallel processes] python3 steady.py -res [resolution] -b [width of cavity] -fx [forcing] -kappa [thermal diffusion coefficient] [...etc]
```
#### Examples
  
To run with the parameters like in Fig. 2a of [Baldelli et al. (2022)], issue the following command:
 
```
mpiexec -np 4 python3 steady.py -res 64 -Lx 120.0 -b 30.0 -H 0.5 -fx 0.01 -Ttop 1e-06 -nutop 0.015 -nubtm 0.015 -kappa 0.000308
```
    
To run with the parameters like in Fig. 2e of [Baldelli et al. (2022)], issue the following command:
 
```
mpiexec -np 4 python3 steady.py -res 64 -Lx 120.0 -b 30.0 -H 0.5 -fx 0.01 -Ttop 8e-06 -nutop 0.015 -nubtm 0.015 -kappa 0.000308
```
  
To run with the parameters like in Fig. 3d of [Baldelli et al. (2022)], issue the following command:
 
```
mpiexec -np 4 python3 steady.py -res 64 -Lx 120.0 -b 30.0 -H 0.5 -fx 0.01 -Ttop 7e-04 -nutop 0.0145 -nubtm 0.0145 -kappa 0.00029
```
  
 ## Visualising the results
The simulations produce the files `U.xdmf` for the velocity field and `T.xdmf` for the temperature field which are best visualized using [ParaView]:
[image of output]

...

The effective heat conduction is computed at the end of each simulation, for example the command
```
(either referenced above or a new one)
```
yields the result `k_eff=X`.
</md>
