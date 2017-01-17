
HOMEDIR = `pwd`



# c compiler
CC 	      = gcc
CC_FLAGS  = -Wall -Wundef -m64 -march=native -O3 -msse4 -ftree-vectorizer-verbose=1 -fopenmp -fPIC
CC_INC    = -I.
CC_DEFS   = 
CC_LIBS   = -lgomp -lpthread -lm -lfftw3

# c++ compiler
CXX	      = g++
CXX_FLAGS = -Wall -Wundef -m64 -march=native -O3 -msse4 -ftree-vectorizer-verbose=1 -fopenmp -fPIC
CXX_INC   = -I.
CXX_DEFS  = -DNX=$(NX) -DNY=$(NY) -DNZ=$(NZ) -DINTERACTIONS=$(INTERACTIONS) -DVEXT=$(VEXT)
CXX_LIBS  = -lgomp -lpthread -lm -lfftw3



FFTW      = fftw_interface.o


all: $(FFTW)

$(FFTW):
	$(CXX) -c fftw_diffs/fftw_interface.c -o $@ $(CXX_FLAGS) $(CXX_INC) $(CXX_DEFS) $(CC_LIBS)
