# vfkm


## Vector Field K-Means Implementation

1) In order to compile, just go to the src directory and make it. This will generate the vfkm binary.

2) The data directory contains sample data files. The data format consists of the first line containing the bouding box of the region used to estimate the vector fields used in the experiment in the form "xmin xmax ymin ymax tmin tmax". Trajectories not contained in this grid are tesselated. The following lines contains the trajectory themselves as a sequence of points and ending with the line "0 0 0".

3) The output of the algorithm is written in the form of tree where the root cluster is represented by "r" and its children are numbered. The vf files contain the vector fields for the cluster and the curve files contain the curves in each cluster.


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