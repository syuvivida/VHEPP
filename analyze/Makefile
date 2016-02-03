# ----------------------------------------
# inputs 

FASTJETLOCATION=/afs/cern.ch/work/s/syu/vlhc/fastjet/fastjet-install/
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

#FILES = Qjets.o QjetsPlugin.o

all: anaSubstructure

anaSubstructure: anaSubstructure.cpp $(FILES)
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lvectorDict -lEnergyCorrelator -lRecursiveTools anaSubstructure.cpp
    
%.o: %.C %.h
	gcc -fPIC -O3 -c $(FASTJETINC) $< -o $@ 

clean:
	rm anaSubstructure