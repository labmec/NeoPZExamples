# NeoPZ Examples
This repository aims to illustrate the functionalities available in the [NeoPZ](https://github.com/labmec/neopz) library.


## Usage
First, install NeoPZ according to the instructions in the GitHub repository.

Then, this project can be configured through CMake. Since some projects rely on graphical output in the `vtk` format, it is recommended to use [Paraview](https://www.paraview.org/) for visualization.

## Examples
### Special Maps

This sample project is used for demonstrating some of the unusual special mappings available in NeoPZ. Relies on graphical output.

### F17 Directional Refinement

 This project used a F17 aircraft as a model to demonstrate the capabilities of automatic directional refinement in NeoPZ. These capabilities might be useful, for instance, in the case of boundary layer analysis, hence the chosen example. The program first builds a mesh modelling the *outside* of the aircraft, and the aircraft is part of the boundary of the created mesh. Then, through sucessive iterations, it models the boundary layer of the jet. Relies on graphical output.