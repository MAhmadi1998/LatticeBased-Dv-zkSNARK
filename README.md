# LatticeBased-Dv-zkSNARK
This is the source code for paper "Lattice-Based Designated Verifier zkSNARKs From Standard
Assumptions"
The Work is based on the source of [ Implementation of lattice trapdoors on modules and appli-
cations]~\ref{https://github.com/lucasprabel/module_gaussian_lattice.
}
To compile the file use the following command
make mainprotocol.c
Then run it with 
./mainprotocol
Note that you can chnage the length of witness by changing the parameter "PARAM_NumTrapdoors" in the "common.h" file
