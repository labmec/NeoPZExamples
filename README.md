# NeoPZ Examples
This repository aims to illustrate the functionalities available in the [NeoPZ](https://github.com/labmec/neopz) library and to provide a short tutorial for introducing the most used classes in the library.


## Usage
First, install NeoPZ according to the instructions in the GitHub repository.

Then, this project can be configured through CMake. Since some projects rely on graphical output in the `vtk` format, it is recommended to use [Paraview](https://www.paraview.org/) for visualization.

## Examples
### Special Maps

This sample project is used for demonstrating some of the unusual special geometrical mappings available in NeoPZ. It also shows how even these curved elements can be geometrically refined. It relies on graphical output.

### F17 Directional Refinement

 This project used a F17 aircraft as a model to demonstrate the capabilities of automatic directional refinement in NeoPZ. These capabilities might be useful, for instance, in the case of boundary layer analysis, hence the chosen example. The program first builds a mesh modelling the *outside* of the aircraft, and the aircraft is part of the boundary of the created mesh. Then, through sucessive iterations, it models the boundary layer of the jet. It relies on graphical output.
 
### ModalAnalysisWG

This project illustrates how to use NeoPZ to perform  modal analysis of electromagnetic Waveguides. It is suitable (as is)for rectangular homogeneous waveguides surrounded with perfect metallic conductors.  It uses as an example a WR90 rectangular waveguide.

It is also an interesting example on how to deal with multiphysics meshes and Krylov eigensolvers with spectral transformations. *This target depends on NeoPZ being linked against LAPACK*.

## Tutorials 
### Poisson2D

Solve the Poisson equation in a bidimensional domain and performs error analysis.

### HCurlProjection

Projects a given analytic solution on a HCurl-conforming approximation space and performs error analysis.