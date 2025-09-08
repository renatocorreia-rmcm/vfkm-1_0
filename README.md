# Vector Field K-Means


1) ### Compilation

    Go to the src directory and run `make`. This will generate the vfkm binary.

2) ### Arguments

    VFKM takes the following: 
    1. `trajectoryFile`
    2. `gridResolution`
    3. `numberOfVectorFields` 
    4. `smoothnessWeight`
    5. `outputDirectory`

    The `data` directory contains sample dataset files ready to use.

3) ### Output
    The output of the algorithm is written in the form of tree where the root cluster is represented by "r" and its children are numbered. The vf files contain the vector fields for the cluster and the curve files contain the curves in each cluster.


## Project Directory tree

```
VFKM-1_0
|   LICENSE
|   README.md
|
+---data                          # DataSets of trajectories ready to use.
|       atlantic_storms.txt
|       msr_tracks.txt
|       msr_tracks_msrcampus.txt
|       synthetic.txt
|       trajectories.txt
|
+---doc
|       instructions.pdf          # More details about this program.
|
\---src                           # Run `make` here.
        ConstraintMatrix.cpp
        ConstraintMatrix.h
        Grid.cpp
        Grid.h
        main.cpp
        Makefile
        Optimizer.cpp
        Optimizer.h
        PolygonalPath.cpp
        PolygonalPath.h
        Util.cpp
        Util.h
        Vector.cpp
        Vector.h
        Vector2D.cc
        Vector2D.h
```