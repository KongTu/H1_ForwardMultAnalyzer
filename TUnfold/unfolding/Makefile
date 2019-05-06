ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTVERSION  := $(shell $(ROOTCONFIG) --version)
ROOTMAJORVERSION := $(word 1,$(subst ., ,$(ROOTVERSION)))

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -fPIC $(ROOTCFLAGS) -I$(SRC)

TUNFOLDDIR=TUnfold_V17.7

LDFLAGS=$(ROOTLDFLAGS) -L. -Wl,-rpath $(TUNFOLDDIR) -Wl,-rpath .
CXX=$(shell $(ROOTCONFIG) --cxx)

ROOTLIBS     := -L$(shell $(ROOTCONFIG) --libdir) -lXMLParser $(ROOTLIBS)

#USE_PRIVATE_TUNFOLD=
#TUNFOLD_LIBS=
USE_PRIVATE_TUNFOLD=-DPRIVATE_TUNFOLD
TUNFOLDVERSION :=$(shell grep TUnfold_VERSION $(TUNFOLDDIR)/TUnfold.h)
TUNFOLDVERSION :=$(filter-out \#define TUnfold_VERSION,$(TUNFOLDVERSION))
TUNFOLDVERSION :=$(shell echo $(TUNFOLDVERSION))
TUNFOLDCODEVER :=$(subst ., ,$(TUNFOLDVERSION))
TUNFOLDCODEVER :=$(word 1,$(TUNFOLDCODEVER))
TUNFOLD_SOURCE=$(shell ls $(TUNFOLDDIR)/TUnfold*.h | grep -v Dict) $(shell ls $(TUNFOLDDIR)/TUnfold*$(TUNFOLDCODEVER).cxx)
TUNFOLD_LIBRARY=$(TUNFOLDDIR)/libunfold.so
TUNFOLD_LIBS=-L$(TUNFOLDDIR) -lunfold

SRC=src
UNFOLDLIBRARY=libmyunfold.so
UNFOLDLIB=-L. -lmyunfold

LIBOBJECTS=Correlations.o UnfoldingResult.o UnfoldingAverage.o UnfoldingAlgorithm.o UnfoldingTUnfoldDensity.o UnfoldingIterative.o VarList.o ClassifierBinning.o
MAKELIBRARY = $(CXX) $(CXXFLAGS) -shared -o $(2) $(1)

all: Makefile fillHist unfoldHist drawResults

$(TUNFOLD_LIBRARY): $(TUNFOLD_SOURCE)
	cd $(TUNFOLDDIR) && make lib

clean:
	rm -f fillHist unfoldHist drawResults unfoldHistmaker $(UNFOLDLIBRARY) \
		fillMerged unfoldMerged drawMerged

%.o: $(SRC)/%.cxx
	$(CXX) $(USE_PRIVATE_TUNFOLD) $(CXXFLAGS) -c $< -o $@

$(UNFOLDLIBRARY): $(LIBOBJECTS)
	$(call MAKELIBRARY,$(LIBOBJECTS),$(@))

fillHist: $(SRC)/fillHist.C $(UNFOLDLIBRARY) $(TUNFOLD_LIBRARY)
	rm -f $@
	$(CXX) $(CXXFLAGS) $(USE_PRIVATE_TUNFOLD) $< -o $@ $(LDFLAGS) \
	$(UNFOLDLIB) $(TUNFOLD_LIBS) $(ROOTLIBS)

unfoldHist: $(SRC)/unfoldHist.C $(UNFOLDLIBRARY) $(TUNFOLD_LIBRARY)
	rm -f $@
	$(CXX) $(CXXFLAGS) $(USE_PRIVATE_TUNFOLD) $< -o $@ $(LDFLAGS) \
	$(UNFOLDLIB) $(TUNFOLD_LIBS) $(ROOTLIBS)

drawResults: $(SRC)/drawResults.C $(UNFOLDLIBRARY) $(TUNFOLD_LIBRARY)
	rm -f $@
	$(CXX) $(CXXFLAGS) $(USE_PRIVATE_TUNFOLD) $< -o $@ $(LDFLAGS) \
	$(UNFOLDLIB) $(TUNFOLD_LIBS) $(ROOTLIBS)
