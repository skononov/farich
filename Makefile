
EXEC := farichres

all: $(EXEC)

CPPFLAGS += -I. -I$(shell root-config --incdir)
CXXFLAGS += -I. $(shell root-config --cflags) $(shell gsl-config --cflags) -fPIC
LD_FLAGS += -Wl,-rpath=$(shell root-config --libdir) -Wl,-rpath-link=$(shell root-config --libdir)
LIBS := $(shell root-config --libs) $(shell gsl-config --libs)
LD := $(CXX)

SOURCES := MLADescription.cc Spectrum.cc farichres.cc
OBJECTS := $(SOURCES:.cc=.o)
DEPS := $(SOURCES:.cc=.cc.d)

%.o: %.cc
	$(CXX) -g -c $(CXXFLAGS) $< -o $@

%.cc.d: %.cc
	@echo "Gererating dependencies for $<"    
	@$(CC) -MM $(CPPFLAGS) $< -MF $@

farichres: $(OBJECTS)
	$(LD) $(LD_FLAGS) $(OBJECTS) $(LIBS) -o $@

ifeq ($(findstring clean,$(MAKECMDGOALS)),)
  -include $(DEPS)
endif

clean:
	rm -f farichres $(OBJECTS)

depclean: 
	rm -f farichres $(OBJECTS) $(DEPS)
