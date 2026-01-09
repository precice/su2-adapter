# SU2-preCICE adapter

This is a [preCICE](https://precice.org/) adapter for [SU2](https://su2code.github.io/) version 7.5.1 "Blackbird",
including examples for conjugate heat transfer and fluid-structure interation simulations.

Read more in the [documentation](https://precice.org/adapter-su2-overview.html).

**Note:** In its current state, the SU2 adapter is using the Python wrapper of SU2. The [previous implementation](https://github.com/precice/su2-adapter/tree/ab843878c1d43302a4f0c66e25dcb364b7787478) was directly implementing the adapter in the C++ source files of SU2. There is also a previous, [unreleased version for preCICE v2](https://github.com/precice/su2-adapter/commit/a87a1ed57e14dca97f1e47aab44632a254714004). If you are looking for the documentation of the previous adapter, get the PDF export of the preCICE website from the [preCICE Distribution v2211.0](https://precice.org/installation-distribution.html#v22110).