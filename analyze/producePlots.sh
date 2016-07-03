python simplePlot_response.py ../data/tev5mumu_mg5_2HDMexovv_rfull006_dat_includeNeutrino 0.4 0 -b
python simplePlot_response.py ../data/tev5mumu_mg5_2HDMexovv_rfull006_dat_includeNeutrino 0.8 0 -b
python simplePlot_response.py ../data/tev5mumu_mg5_2HDMexovv_rfull006_dat_includeNeutrino 1.5 0 -b
mv plots_radius* ../plots/response/2HDM/.

python simplePlot_response.py ../data/tev40mumu_pythia6_zprime5tev_wwrfull006_onlyhadronic 0.4 1 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime5tev_wwrfull006_onlyhadronic 0.8 1 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime5tev_wwrfull006_onlyhadronic 1.5 1 -b
mv plots_radius* ../plots/response/40TeVmumu/zprime5tev/.


python simplePlot_response.py ../data/tev40mumu_pythia6_zprime10tev_wwrfull006_onlyhadronic 0.4 2 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime10tev_wwrfull006_onlyhadronic 0.8 2 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime10tev_wwrfull006_onlyhadronic 1.5 2 -b
mv plots_radius* ../plots/response/40TeVmumu/zprime10tev/.


python simplePlot_response.py ../data/tev40mumu_pythia6_zprime20tev_wwrfull006_onlyhadronic 0.4 3 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime20tev_wwrfull006_onlyhadronic 0.8 3 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime20tev_wwrfull006_onlyhadronic 1.5 3 -b
mv plots_radius* ../plots/response/40TeVmumu/zprime20tev/.


python simplePlot_response.py ../data/tev40mumu_pythia6_zprime30tev_wwrfull006_onlyhadronic 0.4 4 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime30tev_wwrfull006_onlyhadronic 0.8 4 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime30tev_wwrfull006_onlyhadronic 1.5 4 -b
mv plots_radius* ../plots/response/40TeVmumu/zprime30tev/.


python simplePlot_response.py ../data/tev40mumu_pythia6_zprime40tev_wwrfull006_onlyhadronic 0.4 5 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime40tev_wwrfull006_onlyhadronic 0.8 5 -b
python simplePlot_response.py ../data/tev40mumu_pythia6_zprime40tev_wwrfull006_onlyhadronic 1.5 5 -b
mv plots_radius* ../plots/response/40TeVmumu/zprime40tev/.

python simplePlot_response.py ../data/tev10mumu_pythia6_zprime10tev_ww_dat_onlyhadronic 0.4 6 -b
python simplePlot_response.py ../data/tev10mumu_pythia6_zprime10tev_ww_dat_onlyhadronic 0.8 6 -b
python simplePlot_response.py ../data/tev10mumu_pythia6_zprime10tev_ww_dat_onlyhadronic 1.5 6 -b
mv plots_radius* ../plots/response/ZprimeWW/.
