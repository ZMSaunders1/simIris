#include Makefile.arch
#Makefile
BASEDIR = $(shell pwd)
SOURCEDIR = $(BASEDIR)/src
INCLUDEDIR = $(BASEDIR)/include
OBJECTDIR = $(BASEDIR)/obj
BINARYDIR = $(BASEDIR)/bin
HEADER = -I$(INCLUDEDIR) #-I$(TREEIRIS)/include

CXX = g++
LD = g++
ifdef ROOTSYS
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs) -lThread -Wl,-rpath,$(ROOTSYS)/lib
CXXFLAGS += -g -O -Wall -Wuninitialized -I./ -I$(ROOTSYS)/include
ROOTCFLAGS    = $(shell root-config --cflags)
CXXFLAGS += $(HEADER)
CXXFLAGS      += -g -ansi -fPIC $(ROOTCFLAGS)
endif 

SOFLAGS       = -g -shared
LDFLAGS = -O2

all:  $(BINARYDIR)/simIris

$(BINARYDIR)/simIris: $(OBJECTDIR)/simIris.o $(OBJECTDIR)/nucleus.o $(OBJECTDIR)/params.o $(OBJECTDIR)/eloss.o $(OBJECTDIR)/shieldClear.o $(OBJECTDIR)/YYHitDict.o $(OBJECTDIR)/CsIHitDict.o $(OBJECTDIR)/S3HitDict.o
	$(LD) -o $@ $(LDFLAGS) $(ROOTGLIBS) $^ #-lm -lz -lutil 

$(OBJECTDIR)/simIris.o: $(SOURCEDIR)/simIris.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJECTDIR)/params.o: $(SOURCEDIR)/params.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/nucleus.o: $(SOURCEDIR)/nucleus.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/eloss.o: $(SOURCEDIR)/eloss.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/shieldClear.o: $(SOURCEDIR)/shieldClear.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/YYHit.o: $(SOURCEDIR)/YYHit.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/YYHitDict.o: $(SOURCEDIR)/YYHitDict.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOURCEDIR)/YYHitDict.cxx: $(SOURCEDIR)/YYHit.cxx $(INCLUDEDIR)/YYHitLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(HEADER) $^

$(OBJECTDIR)/CsIHit.o: $(SOURCEDIR)/CsIHit.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@
 
$(OBJECTDIR)/CsIHitDict.o: $(SOURCEDIR)/CsIHitDict.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOURCEDIR)/CsIHitDict.cxx: $(SOURCEDIR)/CsIHit.cxx $(INCLUDEDIR)/CsIHitLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(HEADER) $^

$(OBJECTDIR)/S3Hit.o: $(SOURCEDIR)/S3Hit.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTDIR)/S3HitDict.o: $(SOURCEDIR)/S3HitDict.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOURCEDIR)/S3HitDict.cxx: $(SOURCEDIR)/S3Hit.cxx $(INCLUDEDIR)/S3HitLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(HEADER) $^

clean::
	rm -f $(OBJECTDIR)/*.o
	rm -f simIris
	rm -f $(SOURCEDIR)/*Dict.cxx
	rm -f $(SOURCEDIR)/*Dict.h

# end 
