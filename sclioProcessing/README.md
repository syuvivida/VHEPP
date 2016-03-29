```
wget http://atlaswww.hep.anl.gov/asc/jas4pp/download/current.php -O jas4pp.tgz
tar -zvxf jas4pp.tgz
cp mc_singlepi.py jas4pp/.
cd jas4pp
bash
source setup.sh
./fpad mc_singlepi.py

```