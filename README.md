# Modal Analysis

Matlab code for modal analysis of 3D objects.

## Information

```
function [V,F,U,q] = computeModalBasis(filename, varargin)
    % Input:
    %   filename  holding volume or surface mesh (in m)
    %             (.stl, .obj, .msh)
    %             if a volume mesh is provided, the corresponding surface mesh with the same
    %             name should be located next to it
    % Output: 
    %   .msh      holding the generated volume mesh if a surface mesh is provided (generated by tetgen)
    %   .javabin  holding the modal basis (k modes) 
    %
    % Options:
    %   'Young'   followed by Young's modulus (Pa) of object (default 1e4)
    %   'Density' followed by density (kg/m3) of object (default 5000)
    %   'NbModes' followed by number of modes to compute (default 10)
    %   'Indices' followed by list of indices of fixed nodes (default none)
    %   'Boxes'   followed by list of axis aligned boxes ROI for fixed nodes (default none)
    %   'Draw'    followed by true or false to visualize the modes
```

## Requirements
 1- Install tetgen: http://wias-berlin.de/software/tetgen/ 
 
## Examples

`computeModalBasis("models/snake.stl", "Draw", true);`

If the volume mesh already exists just run:

`computeModalBasis("models/snake.msh", "Draw", true);`
