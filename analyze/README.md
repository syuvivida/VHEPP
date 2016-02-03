#Install Fastjet and Fastjet contrib

cd ../
mkdir fastjet
cd fastjet
curl -O http://fastjet.fr/repo/fastjet-3.1.3.tar.gz 
tar zxvf fastjet-3.1.3.tar.gz
cd fastjet-3.1.3/
./configure --prefix=$PWD/../fastjet-install
make 
make check
make install
cd ..
curl -O http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.021.tar.gz
tar xvzf fjcontrib-1.021.tar.gz
cd fjcontrib-1.021
./configure --fastjet-config=/Users/ntran/Documents/Research/Ext/VHEPP/VHEPPStudies/fastjet/fastjet-install/bin/fastjet-config
make 
make check
make install

#Then run the code
cd ../analyze
make
./anaSubstructure of_PanPFA