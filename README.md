# iRED
Scripts to perform iRED (isotropic reorientational eigenmode dynamic) analysis of bond vectors to obtain order parameters (S2). For Amber 14 or later.
iRED analysis is described in Prompers, J. J., Bruchweiler, R. J. Am. Chem. Soc. 2002,124, 4522-4534. DOI: [10.1021/ja012750u](http://dx.doi.org/10.1021/ja012750u)

First,we need to create input to cpptraj that tells which bond vectors that are of interest. This is created by a Python script called `create_bv_inpt.py`. Type,

    python create_bv_inpt.py -v nh -p ref.pdb -o iredout -t mdcrd5 > cpptraj.in

* `nh` - the type of bond vectors of interest. See below, for the supported types of vectors.
* `ref.pdb` - a reference pdb-file used to find the bond vectors. Also used as the topology file.
* `iredout` - a cpptraj output file template string
* `mdcrd5` - the MD trajectory. You can specify one of more, separated by space. If you would like to only process a part of the trajectory, enclose the trajectory and the selection in quotation-marks, e.g., "mdcrd5 1 100"

`cpptraj.in` now contains the following lines (or something similar):

    param ref.pdb 
    trajin mdcrd5
    vector v1 :2@N ired :2@H
    vector v2 :3@N ired :3@H
    ...
    matrix ired name matired order 2
    diagmatrix matired vecs 128 name ired.vec
    ired order 2 modes ired.vec orderparamfile iredout_s2 out iredout

The first two lines is a standard lines that reads in the topology and the MD trajectory. Then follows 128 lines with bond vector definitions, where 128 is the total number of bond vectors. In the example above,we define backbone N-H bond vectors. The last three lines perform the actual creation of the covariance matrix, the calculation of the eigenvalues/eigenvectors, and the calculation of the order parameters, respectively. Here it is important that vecs are set to the total number of bond vectors.

Now, run cpptraj

    cpptraj prmtop cpptraj.in >& cpptraj.out

