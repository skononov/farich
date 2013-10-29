
EXEC := farichres

all: $(EXEC)

CXXFLAGS += -I. $(shell root-config --cflags) $(shell gsl-config --cflags)
LD_FLAGS += -lMinuit2 $(shell root-config --libs) $(shell gsl-config --libs)
LD := $(CXX)

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $^ -o $@

farichres: farichres.o MLADescription.o Spectrum.o Minuit2FCN.o
	$(LD) $^ $(LD_FLAGS) -o $@

clean:
	rm -f farichres *.o
