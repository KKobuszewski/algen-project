
HOMEDIR = `pwd`

GPUARCH = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52

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


all: eigensolver hamiltonian swap simmang harmonic
#$(FFTW)

$(FFTW):
	$(CXX) -c fftw_diffs/fftw_interface.c -o $@ $(CXX_FLAGS) $(CXX_INC) $(CXX_DEFS) $(CC_LIBS)

eigensolver:
	nvcc -c $(HOMEDIR)/eigensolver/Eigensolver.cu -o $(HOMEDIR)/eigensolver/Eigensolver.o -O3 $(GPUARCH) -Xcompiler "-fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm"
	nvcc -c $(HOMEDIR)/tests/test_eigensolver.o $(HOMEDIR)/tests/test_eigensolver.cpp $(GPUARCH) -Xcompiler "-std=c++11 -fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm"
	nvcc -o $(HOMEDIR)/tests/test_eigensolver.exe $(HOMEDIR)/tests/test_eigensolver.o $(HOMEDIR)/eigensolver/Eigensolver.o $(GPUARCH) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -L/usr/local/cuda/lib64 -lcusolver -lopenblas -lgomp -lm -lcudart

hamiltonian:
	@echo $(HOMEDIR)
	g++ -std=c++11 $(HOMEDIR)/tests/test_hamiltonian.cpp -o $(HOMEDIR)/tests/test_hamiltonian.exe -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -lgomp -lm -DHARMONIC_POTENTIAL

swap:
	nvcc -c $(HOMEDIR)/eigensolver/Eigensolver.cu -o $(HOMEDIR)/eigensolver/Eigensolver.o -O3 $(GPUARCH) -Xcompiler "-fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm"
	g++ -std=c++11 -c $(HOMEDIR)/tests/test_swap_basis.cpp -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -I/usr/local/cuda/include -DHARMONIC_POTENTIAL -DTOURNAMENT_SELECTION
	nvcc -o $(HOMEDIR)/tests/test_swap_basis.exe test_swap_basis.o $(HOMEDIR)/eigensolver/Eigensolver.o $(GPUARCH) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -L/usr/local/cuda/lib64 -lcusolver -lopenblas -lgomp -lm -lcudart
	@rm test_swap_basis.o

simmang:
	nvcc -c $(HOMEDIR)/eigensolver/Eigensolver.cu -o $(HOMEDIR)/eigensolver/Eigensolver.o -O3 $(GPUARCH) -Xcompiler "-fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm"
	g++ -std=c++11 -c $(HOMEDIR)/tests/test_simmang.cpp -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -I/usr/local/cuda/include -DHARMONIC_POTENTIAL -DTHREADS=5
	nvcc -o $(HOMEDIR)/tests/test_simmang.exe test_simmang.o $(HOMEDIR)/eigensolver/Eigensolver.o $(GPUARCH) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -L/usr/local/cuda/lib64 -lcusolver -lopenblas -lgomp -lm -lcudart
	@rm test_simmang.o

harmonic:
	nvcc -c $(HOMEDIR)/eigensolver/Eigensolver.cu -o $(HOMEDIR)/eigensolver/Eigensolver.o -O3 $(GPUARCH) -Xcompiler "-fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm"
	g++ -std=c++11 -c $(HOMEDIR)/tests/test_harmonic.cpp -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -I/usr/local/cuda/include -DHARMONIC_POTENTIAL -DRANDOM_SELECTION
	nvcc -o $(HOMEDIR)/tests/test_harmonic_rnd_select.exe test_harmonic.o $(HOMEDIR)/eigensolver/Eigensolver.o $(GPUARCH) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -L/usr/local/cuda/lib64 -lcusolver -lopenblas -lgomp -lm -lcudart
	@rm test_harmonic.o
	g++ -std=c++11 -c $(HOMEDIR)/tests/test_harmonic.cpp -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -I/usr/local/cuda/include -DHARMONIC_POTENTIAL -DTOURNAMENT_SELECTION
	nvcc -o $(HOMEDIR)/tests/test_harmonic_tournament_select.exe test_harmonic.o $(HOMEDIR)/eigensolver/Eigensolver.o $(GPUARCH) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -L/usr/local/cuda/lib64 -lcusolver -lopenblas -lgomp -lm -lcudart
	@rm test_harmonic.o

debug:
	g++ -std=c++11 -g -Og -ggdb $(HOMEDIR)/tests/test_simmang.cpp -o $(HOMEDIR)/tests/test_simmang_debug.exe -fPIC -fopenmp -mtune=native -march=native -O3 -I$(HOMEDIR) -L$(HOMEDIR)/tests -L$(HOMEDIR)/eigensolver -lopenblas -lgomp -lm -lc -DHARMONIC_POTENTIAL -DTHREADS=5
	#g++ -std=c++11 -g -Og -ggdb test_hamiltonian.cpp -o test_hamiltonian_debug.exe -fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm -DHARMONIC_POTENTIAL

clean:
	@rm *.o
	@rm $(HOMEDIR)/eigensolver*.o
	@rm *.exe