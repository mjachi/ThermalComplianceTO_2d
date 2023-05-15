# Topology Optimization

MFEM topology optimization code for thermal compliance.

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

