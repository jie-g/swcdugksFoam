## About

This project has successfully ported the `cdugksFoam` solver to the Sunway TaihuLight platform and made a series of targeted optimizations.

## Installation

Before installing the code in this repository, we assume that you have installed the MPI library and swOpenFOAM.

Then, installation is carried out according to the following procedure:

```
git clone https://github.com/jie-g/swcdugksFoam.git
cd swcdugksFoam/src/fvDVM/swFunctions/src
./run
cd swcdugksFoam/src/
./Allwmake
```

Tip: You can change the macros in `para.h` to suit your needs. It will produce a different effect. Note that you need to recompile after changing.
## Demo

The 2D lid-driven cavity flow simulation in the paper is provided in `swcdugksFoam/demo/cavity0.1`, which you can try to test.
## References

- [1] L. Zhu, [dugksFoam](https://github.com/zhulianhua/dugksFoam).  
- [2] Zhang Q, [cdugksFoam](https://github.com/zzhang777/paralleled_cdugksFoam/tree/master).  
- [3] Zhang Q, Wang Y, Pan D, Unified X-space parallelization algorithm for conserved discrete unified gas kinetic scheme. [Computer Physics Communications, 2022, 278: 108410](https://www.sciencedirect.com/science/article/abs/pii/S0010465522001291).  
- [4] L. Zhu, S. Chen, Z. Guo, dugksFoam: An open source OpenFOAM solver for the Boltzmann model equation, [Comp. Phys. Commun., 213(2017) 155-164](http://www.sciencedirect.com/science/article/pii/S0010465516303642).  
- [5] Z. Guo, K. Xu, R. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305).  
- [6] Z. Guo, R. Wang, K. Xu, Discrete unified gas kinetic scheme for all Knudsen number flows. II. Thermal compressible case, [Phys. Rev. E, 91(2015) 033313](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.033313).  
- [7] L. Zhu, Z. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, [Comp. Fluids, 127(2016) 211-225](http://www.sciencedirect.com/science/article/pii/S0045793016000177).  
- [8] Fu H, Liao J, Yang J,The Sunway TaihuLight supercomputer: system and applications, [Science China Information Sciences, 2016, 59: 1-16](https://link.springer.com/article/10.1007/s11432-016-5588-7).  
