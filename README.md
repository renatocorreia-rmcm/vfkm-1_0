# Vector Field K-Means

1) ### Compilation
        Go to the src directory and run `make`. This will generate the vfkm binary.

2) ### Data
        The `data` directory contains sample dataset files. See the description of its structure in `doc/instructions.pdf` - `section 3`.

3) ### Output
        The output of the algorithm is written in the form of tree where the root cluster is represented by "r" and its children are numbered. The vf files contain the vector fields for the cluster and the curve files contain the curves in each cluster.


## Directory tree

```
C:.
|   LICENSE
|   README.md
|
+---data
|       atlantic_storms.txt
|       msr_tracks.txt
|       msr_tracks_msrcampus.txt
|       synthetic.txt
|       trajectories.txt
|
+---doc
|   |   instructions.pdf
|   |   instructions.tex
|   |
|   \---figs
|           file_example.pdf
|           first_example.pdf
|           grid_example.pdf
|           vertex_ordering.pdf
|
\---src
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