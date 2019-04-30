
CXXFLAGS := $(shell root-config --cflags) -I$(HOME)/root/include -I. -O2
LDOPTS := -O $(shell root-config --ldflags) -L$(HOME)/root/lib -lMyRoot
LIBS := $(shell root-config --libs)

SOURCES := ctaulayerringcal.cc

.PHONY: all clean

all: richcal

%.o: %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

richcal: $(SOURCES:.cc=.o)
	$(CXX) $(LDOPTS) -o $@ $^ $(LIBS)

clean:
	rm -f $(SOURCES:.cc=.o) richcal
