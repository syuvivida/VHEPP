## Download codes from github

## run the code
```
cd ../analyze
make -f Make.rawhit
./anaRawhit dataDirectory 0.4 0 0.5 // for mode 0, cut E>0.5 GeV
./anaRawhit dataDirectory 0.4 1 90 // for mode 1, remove the lowest 90% of hits
```