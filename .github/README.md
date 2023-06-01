# Topology Optimization: Thermal Compliance

MFEM topology optimization code for thermal compliance, implementing 5 simple problems.

### Building

The serial code can be built with CMake as follows. If MFEM is already built somewhere that it can find, it should recognize it;
otherwise, build MFEM and GLVis in e.g. a `lib` folder under the root as usual.
```
> mkdir build && cd build
> cmake ..
> make
> ./main
```
The CMake does not build the parallel example, however, `mainp.cxx` and `main.hxx` can be copy/pasted into the `examples` folder
within the `mfem-*.*.*`, etc, directory and added to the parallel examples build rules. The file extension should be changed to
`mainp.cpp` for compatibility.

#### Warning

For several of the problems, the default values for certain parameters may not allow for convergence of the method. In particular,
if you see `Projection reached maximum iteration without converging. Result may not be accurate.`, you might consider lowering
the step size (the argument `-alpha`) e.g. by a factor of 10/ order of magnitude. Lowering the mass fraction (the CLI argument `-mf`)
to `0.3` or `0.4` can assist this and also lead to more qualitatively interesting results, particularly on Problem 5.
