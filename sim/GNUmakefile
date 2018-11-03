name:=aRichSim

G4TARGET = $(name)
G4EXLIB := 1

G4UI_USE_TCSH := 1
G4DEBUG := 1
G4LIB_NO_SHARED := 1

DEFINES:=

include $(G4INSTALL)/config/architecture.gmk

.PHONY: default ms bg bge full libroot

default: lib bin

CPPFLAGS += `root-config --cflags` -I$(HOME)/root/include -DWORKDIR=\"$(PWD)\" $(DEFINES)
CXXFLAGS += -pthread -Wno-unused -O3

EXTRALIBS := -L$(RHOME)/lib -Xlinker -rpath -Xlinker $(RHOME)/lib -lMyRoot
EXTRALIBS += -Xlinker -rpath -Xlinker $(G4LIB)/$(G4SYSTEM)
EXTRALIBS += -Xlinker -rpath -Xlinker $(CLHEP_LIB_DIR)
EXTRALIBS += `root-config --libs` -Xlinker -rpath -Xlinker `root-config --libdir`

#CPPVERBOSE := 1

include $(G4INSTALL)/config/binmake.gmk

#all: default ms bg bge full
all: default ms bg full libroot

RichRootDict.h RichRootDict.cc: include/RichRootEvent.hh include/RichPrimaryData.hh include/RichLinkDef.hh
	@echo "Generating ROOT dictionary for RichPrimaryData"
	@rootcint -f $@ -c $^

$(G4TMPDIR)/RichRootDict.o: RichRootDict.cc
	@echo Compiling RichRootDict.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(OUT_OBJ)$@ $<

$(G4LIBDIR)/lib$(name).a: $(G4TMPDIR)/RichRootDict.o

libroot: libroot.so

libroot.so: $(G4TMPDIR)/RichPrimaryData.o $(G4TMPDIR)/RichRootEvent.o $(G4TMPDIR)/RichRootDict.o
	@echo Building shared library for cint...
	@$(CXX) -shared $^ -o $@

ms:
	@echo "Making simulation with MS on"
	@$(MAKE) name=aRichSimMS DEFINES='-DMS_ON=1'

bg:
	@echo "Making simulation with MS and Background on"
	@$(MAKE) name=aRichSimBg DEFINES='-DMS_ON=1 -DBACKGROUND_ON=1'

bge:
	@echo "Making simulation of electron background"
#	@$(MAKE) name=aRichSimBgE DEFINES='-DBGE_ON=1 -DBACKGROUND_ON=1 -DSECONDARIES_ON=1'
	@$(MAKE) name=aRichSimBgE DEFINES='-DBGE_ON=1 -DMS_ON=1 -DBACKGROUND_ON=1 -DSECONDARIES_ON=1'

full:
	@echo "Making full simulation (MS,Background,Decays - on)"
	@$(MAKE) name=aRichSimAll DEFINES='-DMS_ON=1 -DBACKGROUND_ON=1 -DSECONDARIES_ON=1'

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_* *.heprep

