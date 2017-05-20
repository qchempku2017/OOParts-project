#/bin/sh
g++ main.cpp scflib.cpp intlib.cpp diislib.cpp -o OOPart_RHF -L/home/izzy/opts/libcint/lib64 -lcint -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl -m64 -I${MKLROOT}/include
