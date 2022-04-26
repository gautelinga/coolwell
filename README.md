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
To run with the parameters like in Fig. 3d of [Baldelli et al. (2022)], issue the following command:

```
mpiexec -np 4 python3 steady.py -res 64 -Lx 120.0 -b 30.0 -H 0.5 -fx 0.01 -Ttop 7e-04 -nutop 0.0145 -nubtm 0.0145 -kappa 0.00029
```  
  
To run with the parameters like in the [image](#example_image) below, issue the following command:
 
 <a name="example_code">
    
```
mpiexec -np 4 python3 steady.py -res 64 -Lx 120.0 -b 30.0 -H 0.5 -fx 0.001 -Ttop 1e-05 -nutop 0.0015 -nubtm 0.0015 -kappa 0.0015
```
</a>
  
 ## Visualising the results
The simulations produce the files `U.xdmf` for the velocity field and `T.xdmf` for the temperature field which are best visualized using [ParaView](https://www.paraview.org/). The image below is created in Paraview from the `U.xdmf` and `T.xdmf` files that result from running [this code](#example_code). It shows the temperature field and the stream tracer filter applied to the velocity field.
  
 <a name="example_image">
          
<img src="https://github.com/gautelinga/coolwell/blob/96aecec74d5329c62e82e57a6f4d502d5d3771aa/Images/Tfield_streamlines_RT0.001.png" width="300">
                                                                                                  </a>
  
...

The code also produces the file `parameters.dat` which contains a few quantities calculated from the temperature and velocity fields:
  | quantity | name in `parameters.dat` | description | name in article [^1]|
  |----------|--------------------------|-------------|----------------------|
  | &lt;&#124;<strong>u</strong>&#124;&gt;<sub>well | `u_b` | magnitude of the velocity averaged over the square well | U |
  | &lt;&#124;<strong>u</strong>&#124;&gt;<sub>chnl |`u_0`| magnitude of the velocity averaged over the channel | n/a |
  | Re<sub>well</sub> | `Re_b` | Reynold number in the well | Re |
  | Re<sub>chnl</sub> | `Re_0` | Reynold number in the channel | n/a |
  | Ra | `Ra` | Rayleigh number | n/a |
  | Pe<sub>well | `Pe` | Péclet number in the well | Pe |
  | R<sub>T | `RT` | ratio between driving and buoyancy force | n/a |
  | Nu | `Nu` | Nusselt number at the bottom surface | Nu |
  | Q | `qy` | total heat flux density across the bottom surface | n/a |
  
  
  The file `parameters.dat` also contains the full line of code needed to run the simulation that created it.
  
</md>
