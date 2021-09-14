.DELETE_ON_ERROR:

ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS   := $(shell $(ROOTCONFIG) --glibs)

CXX         := g++
CXXFLAGS    := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD          := g++
LDFLAGS     := -O2 $(ROOTLDFLAGS) 

LIBS        := $(ROOTGLIBS)
DEFINES     := -DWITHRINT

SRC_FILES   := $(wildcard *.cxx)
EXE_FILES   := $(SRC_FILES:%.cxx=%) 

.PHONY: all clean

all: $(EXE_FILES)


%: %.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@


clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)
