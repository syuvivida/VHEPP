#Install Fastjet and Fastjet contrib

## For FastJet
```
cd analyze/
mkdir fastjet
cd fastjet
curl -O http://fastjet.fr/repo/fastjet-3.1.3.tar.gz 
tar zxvf fastjet-3.1.3.tar.gz
cd fastjet-3.1.3/
./configure --prefix=$PWD/../fastjet-install
make && make check && make install
```

##For FastJet/contrib:
```
cd ..
svn checkout http://fastjet.hepforge.org/svn/contrib/trunk fjcontrib
cd fjcontrib/
scripts/update-contribs.sh 
scripts/update-contribs.sh EnergyCorrelator 1.2.0-rc1
./configure --fastjet-config=$PWD/../fastjet-install/bin/fastjet-config
make && make check && make install
```
#Then run the code
```
cd ../analyze
make
./anaSubstructure of_PanPFA
```