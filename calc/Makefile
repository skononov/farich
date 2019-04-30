SOURCEDIR := $(dir $(lastword ${MAKEFILE_LIST}))

EXEC := farichres

all: $(EXEC)

CXXFLAGS += -I${SOURCEDIR} $(shell root-config --cflags) $(shell gsl-config --cflags) -fPIC 
LD_FLAGS += -Wl,-rpath=$(shell root-config --libdir) -Wl,-rpath-link=$(shell root-config --libdir)
LIBS := $(shell root-config --libs) $(shell gsl-config --libs) -lgomp
LD := $(CXX)

SOURCES := MLADescription.cc Spectrum.cc farichres.cc
OBJECTS := $(SOURCES:.cc=.o)
DEPS := $(SOURCES:.cc=.cc.d)

ifdef DEBUG
  CXXFLAGS += -g
else
  CXXFLAGS += -O3 -fopenmp
endif

%.o: ${SOURCEDIR}%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.cc.d: ${SOURCEDIR}%.cc
	@echo "Generating dependencies for $<"    
	@$(CC) -MM $(CXXFLAGS) $< -MF $@

farichres: $(OBJECTS)
	$(LD) $(LD_FLAGS) $(OBJECTS) $(LIBS) -o $@

ifeq ($(findstring clean,$(MAKECMDGOALS)),)
  -include $(DEPS)
endif

clean:
	rm -f farichres $(OBJECTS)

depclean: 
	rm -f farichres $(OBJECTS) $(DEPS)
