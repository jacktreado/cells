# Github repository for `cells` 

Code base for 2D and (in the future) 3D simulations of deformable particles, or cells

written in `c++`

## `src` , `main` and `bash`

Object-oriented framework, slower but contains `hopper` and `jamming` code in `src/` directory

Header files forthcoming.

## sequential

Sequential `c++` code, all single `main` file with functions defined in the `.cpp` file.

For the most part, you can compile with `g++`:

`g++ -Wall -O3 sequential/[CODE_NAME_HERE].cpp -o [OUT_NAME_HERE]` 

### Code that needs additional headers
* `bidRepulsiveCellJamming.cpp`
* `bidRepulsiveCellVDOS.cpp`

These functions require use of the `Eigen` header to compute eigenvalues of the dynamical matrix.

To compile, use:

`g++ -Wall -O3 -I src sequential/[JAMMING_CODE_NAME_HERE].cpp -o [OUT_NAME_HERE]`

