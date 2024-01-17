# SurfgenBound 2 PES analysis tools
This is source code that generates various executables that can be used to analyzes surfaces generated by [SurfgenBound 2](https://github.com/cavanes1/SurfgenBound2). The initial version was provided by Yafu Guan.

The provided makefile runs successfully on ARCH.

## Utilities

### Required input files for all
* basis.in
* Hd.CheckPoint
* refgeom
* intcfl

### dat.x
Required input files
* geom.all
* energy.all
* fit.in
* names.all

### findcp.x
Required input files
* findcp.in
* Geometry file (ex. geom)

For example, perform search on state 1 using "geom" as starting guess by running `findcp.x geom 1`

### findmex.x
Required input files
* mexopt.in
* Geometry file (ex. geom)

For example, perform search on states 1 and 2 using "geom" as starting guess by running `findmex.x geom 1 1`

### gf.x
Required input files
* energy.all
* geom.all

Output file: hess

### basis.x
Required input files
* basisest.in - has one parameter: NVS_contour
* precursor.xyz
* precursor.hess
* residual.xyz
* residual.hess

Output files:
* InternalCoverage.txt
* NormalCoverage.txt
* basis.txt

### nadvibs.x
Required input files
* precursor.xyz
* precursor.hess
* residual.xyz
* residual.hess

Output file: nadvibs.in

Note: Precursor is the anion, residual is what's left after electron is removed

### traj.x
Need more information on what this does.

Required input files
* traj.all