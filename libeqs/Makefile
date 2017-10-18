CC  = gcc
F90 = gfortran
PYTHON = python
OPENMP = -fopenmp
LIB_PATH = ../lib

NAME = eqs

ifeq ($(OS),Windows_NT)
	DLL_PATH = ../x64/Release
	BLASLAPACK = $(DLL_PATH)/libblas.dll $(DLL_PATH)/liblapack.dll
	# BLASLAPACK=libopenblas.dll
	FFTW = $(DLL_PATH)/libfftw3-3.dll
	OBJ_EXT=obj
	DLL_EXT=dll
	LIB_EXT=lib
	DEP_LIB=$(NAME).$(LIB_EXT)
	F90FLAG = -shared -Ofast $(OPENMP)
	PREFIX = ../x64/Release
else
	BLASLAPACK = -lblas -llapack
	FFTW = -lfftw3
	OBJ_EXT =o
	LIB_EXT=a
	DEP_LIB =
	PREFIX = /usr/lib
	ifeq ($(shell uname),Darwin)
		DLL_EXT=dylib
		DLLFLAG=-dynamiclib
		F90FLAG= -fPIC -shared -Ofast
	else
		DLL_EXT=so
		DLLFLAG=
		F90FLAG= -fPIC -shared -Ofast
	endif
endif

DYLIB = lib$(NAME).$(DLL_EXT)
STLIB = lib$(NAME).$(LIB_EXT)
OBJ = $(NAME).$(OBJ_EXT)

all: lib

lib: $(DYLIB) $(DEP_LIB)

dylib: $(DYLIB)

stlib: $(STLIB)

$(DYLIB): $(OBJ)
	$(F90) $(OBJ) $(F90FLAG) $(DLLFLAG) $(BLASLAPACK) $(FFTW) -o $(DYLIB)

$(STLIB): $(OBJ)
	$(AR) rcs $(STLIB) $(OBJ)

%.$(OBJ_EXT):%.f90
	$(F90) -c $< $(F90FLAG) -o $@

%.lib:lib%.dll
	$(PYTHON) dll2lib.py $<

install: $(DYLIB) 
	cp $(DYLIB) $(PREFIX)
ifeq ($(OS),Windows_NT)
	cp $(DEP_LIB) ../EQSignal/lib
else
	cp $(DYLIB) ../EQSignal/lib
endif

uninstall:
	rm -vf $(PREFIX)/$(NAME).$(DLL_EXT)

.PHONY : clean, test

clean:
	rm -rvf *.so *.dll *.dylib *.o *.obj *.lib *.def *.mod __pycache__

test: eqspy.py $(DYLIB)
	python eqspy.py
