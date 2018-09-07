
EXEC := farichres

all: $(EXEC)

CXXFLAGS += -I. $(shell root-config --cflags) $(shell gsl-config --cflags) -fPIC
LD_FLAGS += -Wl,-rpath=$(shell root-config --libdir) -Wl,-rpath-link=$(shell root-config --libdir)
LIBS := $(shell root-config --libs) $(shell gsl-config --libs)
LD := $(CXX)

%.o: %.cc
	$(CXX) -g -c $(CXXFLAGS) $^ -o $@

farichres: farichres.o MLADescription.o Spectrum.o
	$(LD) $(LD_FLAGS) $^ $(LIBS) -o $@

clean:
	rm -f farichres *.o
