CXX=g++
ODIR=obj
BDIR=bin
SDIR=source
IDIR=include

BPP_INCLUDE=/include
BPP_LIB=/lib

BOOST_INCLUDE=/usr/include
BOOST_LIB=/usr/lib



CPPFLAGS= -O3 -I$(IDIR) -I$(BPP_INCLUDE) -I$(BOOST_INCLUDE) -std=c++11
LDFLAGS=-L$(BPP_LIB) -L$(BOOST_LIB)
LIBS=-lm -lbpp-core -lbpp-seq -lbpp-phyl


DECO_FILES = DeCoSTAR.cpp AdjMatrix_CostFunctions.cpp CladesAndTripartitions.cpp DeCoUtils.cpp EquivalenceClassFamily.cpp MyGeneTree.cpp MyMatrixV.cpp XMLUtils.cpp AdjMatrix.cpp CoEvent.cpp DTLGraph.cpp GeneFamily.cpp MyMatrixAdjSolution.cpp MySpeciesTree.cpp AdjTree.cpp DeCoOutputManager.cpp DTLMatrix.cpp MultiRootEquivalenceClass.cpp MyMatrix.cpp ReconciledTree.cpp CladeReconciliation.cpp EquivalenceClass.cpp MyCladesAndTripartitions.cpp MyMatrixT.cpp ReconciliationEvent.cpp
DECO_OBJS = $(subst .cpp,.o,$(DECO_FILES))
DECO_OBJS := $(addprefix $(ODIR)/, $(DECO_OBJS))
DECO_SRCS = $(addprefix $(SDIR)/, $(DECO_FILES))

_EXES= DeCoSTAR
EXES=$(patsubst %,$(BDIR)/%,$(_EXES))


all: $(EXES)

debug: CPPFLAGS = -O3 -I$(IDIR) -I$(BPP_INCLUDE) -I$(BOOST_INCLUDE) -g -std=c++11
debug: all


$(BDIR)/DeCoSTAR: $(DECO_OBJS)


#general rules

$(ODIR)/%.o: $(SDIR)/%.cpp
	mkdir -p $(ODIR)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -c -MD -o $@ $<
	@cp $(ODIR)/$*.d $(ODIR)/$*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $(ODIR)/$*.d >> $(ODIR)/$*.P; \
	rm -f $(ODIR)/$*.d

$(BDIR)/%: $(ODIR)/%.o
	mkdir -p bin 
	$(CXX) -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LIBS)



all: $(EXES)


clean:
		rm -f $(ODIR)/*.o *~ $(SDIR)/*~ core  
		rm -rf $(ODIR) 
