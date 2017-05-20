 # OOParts-project
OOParts: Oops! I hope you know what you are doing!
Version 1.0.0
This version of OOParts still needs to specify the atom and the basis information in data/atom.txt, bas.txt, env.txt manually.
Only support close-shell restricted hatree fock. Now enables DIIS.

Group leader: Fengyu Xie(Izzy)
Acknowledgements:

Qiming Sun: The author of libcint, who had given me so many helpful suggestions during this course project!
Chao Huang: For his excellent humour as a TA and his benevolent tolerance to our procrastination!
Yu Jin & Shiqi Chen: The author of diis library.
Jinsong Zhou: Co-author of scf library.
Hexu Liu & Chaolun Sun: Compiled a simple 1s-GTO library, and implemented the project with MATLAB. 
Wenjian Liu: The mentor of quantum chemistry class, 2017, PKUCCME.

installation:

1,Make sure you have installed intel MKL library;
2,Uncompress and install Libcint correctly. Rememeber adding libcint library path to your .bashrc;
3,Edit compile.sh, change the paths inside in correspondence with your own computer.
4,./compile.sh
5,Set up information in bas.txt,atom.txt,env.txt,initial_c.txt correctly according to libcint documentation. 
6,Run OOParts. If you want to do the calculation again in the same directory, run cleanse.sh and restore the calculation to start point. If you want to continue an unfinished calculation, run OOParts directly.
