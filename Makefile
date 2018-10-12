
EXEC := farichres

all: $(EXEC)

CXXFLAGS += -I. $(shell root-config --cflags) $(shell gsl-config --cflags) -fPIC -O3 #-g
LD_FLAGS += -Wl,-rpath=$(shell root-config --libdir) -Wl,-rpath-link=$(shell root-config --libdir)
LIBS := $(shell root-config --libs) $(shell gsl-config --libs)
LD := $(CXX)

SOURCES := MLADescription.cc Spectrum.cc farichres.cc
OBJECTS := $(SOURCES:.cc=.o)
DEPS := $(SOURCES:.cc=.cc.d)

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.cc.d: %.cc
	@echo "Gererating dependencies for $<"    
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
