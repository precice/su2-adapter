---
title: The SU2 adapter
permalink: adapter-su2-overview.html
redirect_from: adapter-su2.html
keywords: adapter, su2, development, modules
summary: "Modify native SU2 files to couple with other solvers or SU2 itself"
---

The [SU2 adapter](https://github.com/precice/su2-adapter) leverages the SU2 Python wrapper and preCICE Python bindings to couple SU2 using preCICE in a minimally invasive way. The adapter simply updates existing functions and implements new ones in the SU2 Python wrapper that allow for simple preCICE implementation with implicit coupling. The adapter currently works for SU2 version 7.5.1 "Blackbird".

The adapter can simulate the flow domain in both conjugate heat-transfer and fluid-structure interaction applications. It supports:

- Temperature (read/write)
- Heat flux (read/write)
- Force (write)
- Displacement (read)

while the Python scripts provided in the `run/` directory can easily be adapted for more fields.

This adapter has been designed to work when using the compressible solver for unsteady problems with dual-time stepping, for single-zone problems. Implicit coupling currently saves the flow solution, turbulence solution, and the mesh solution (for mesh deformation). Species transport and transition model variables at this time are not saved, but may be straightforward to implement.

{% note %}
In its current state, the SU2 adapter is using the Python wrapper of SU2. The [previous implementation](https://github.com/precice/su2-adapter/tree/ab843878c1d43302a4f0c66e25dcb364b7787478) was directly editing the C++ source files of SU2. There is also a [version relying on the Python wrapper that however works with preCICE v2](https://github.com/precice/su2-adapter/commit/a87a1ed57e14dca97f1e47aab44632a254714004). If you are looking for the documentation of the previous adapter, get the PDF export of the preCICE website from the [preCICE Distribution v2211.0](https://precice.org/installation-distribution.html#v22110).
{% endnote %}

## Try

Here you will find how to [get the adapter](adapter-su2-get.html) and [how to configure](adapter-su2-configure.html) a case.
Tutorial cases that use this adapter include the [perpendicular flap](tutorials-perpendicular-flap.html) and the [flow over a heated plate](tutorials-flow-over-heated-plate.html).

## References

[1] Alexander Rusch. [Extending SU2 to Fluid-Structure Interaction via preCICE](http://www5.in.tum.de/pub/Rusch2016_BA.pdf). Bachelor's thesis, Munich School of Engineering, Technical University of Munich, 2016.
