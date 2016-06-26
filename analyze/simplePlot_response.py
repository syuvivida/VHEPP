
import ROOT
import tdrstyle
import math
import os

import sys
for arg in sys.argv: 
    print arg

inputFolder = sys.argv[1]
jetRadius = sys.argv[2]


tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gStyle.SetPalette(1)

def makeCanvas(hists, tags):
    colors = [1,2,4,8,6,7,8]

#    if hists[0].GetName() == "h_jp_gen_nonu" or hists[0].GetName() == "h_jpt_gen_nonu" or hists[0].GetName() == "h_je_gen_nonu":            
#        leg = ROOT.TLegend(0.4,0.6,0.8,0.9)
#    else:
    if hists[0].GetName() == "h_jeratio_gen_nonu":
        leg = ROOT.TLegend(0.18,0.55,0.8,0.9)
    else:
        leg = ROOT.TLegend(0.18,0.6,0.58,0.9)
    leg.SetHeader("Jet radius="+jetRadius)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    tmax = -999


    for i in range(len(hists)):
        if hists[0].GetName() == "h_jeratio_gen_nonu":
            tagname=tags[i] + ", RMS = " + "{:.3f}".format(hists[i].GetRMS())
            leg.AddEntry(hists[i],tagname,'l')
        else:
            leg.AddEntry(hists[i],tags[i],'l')
        if hists[i].GetMaximum() > tmax:
            tmax = hists[i].GetMaximum()

    if hists[0].GetName() == "h_jp_gen_nonu" or hists[0].GetName() == "h_jpt_gen_nonu" or hists[0].GetName() == "h_je_gen_nonu":
        hists[0].SetMaximum(tmax*1.2)
    else:
        hists[0].SetMaximum(tmax*2.0)
            
    c = ROOT.TCanvas("c","c",1000,800)
    for i in range(len(hists)):
        hists[i].SetMinimum(0)
        hists[i].GetYaxis().SetDecimals()
        hists[i].SetLineWidth(3)
        hists[i].SetLineStyle(i+1)
        if i == 0:
            hists[i].DrawNormalized('hist')
        else: 
            hists[i].SetLineColor( colors[i] )
            hists[i].DrawNormalized('histsames')

        leg.Draw()



        c.SaveAs(directory+"/radius"+jetRadius+"_"+hists[0].GetName()+".pdf")
        c.SaveAs(directory+"/radius"+jetRadius+"_"+hists[0].GetName()+".png")

## tin, PF, gen
def getHists(tin,tinmc,postfix):
## Z' = 40 TeV
#	h_je       = ROOT.TH1F("h_je"+postfix,"; jet energy; N",100,1000,6000);
#	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,6000)
#	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,1000,6000)

        h_jeratio  = ROOT.TH1F("h_jeratio"+postfix,"; E_{jet}/E_{jet}^{true}; N",50,0.5,1.2);
#	h_je       = ROOT.TH1F("h_je"+postfix,"; jet energy; N",50,0,22000);
#	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",50,0,22000)
#	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",50,0,22000)

#	h_je       = ROOT.TH1F("h_je"+postfix,"; jet energy; N",50,-15000,22000);
#	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",50,0,22000)
#	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",50,-15000,22000)

# 2HDM
#       h_jeratio  = ROOT.TH1F("h_jeratio"+postfix,"; jet energy/0.5 M_{resonance}; N",50,0.2,1.2);
	h_je       = ROOT.TH1F("h_je"+postfix,"; jet energy; N",100,1000,3000);
	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,3000)
	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,1000,3000)

        h_jeratio.SetNdivisions(5);
        h_je.SetNdivisions(5);
        h_jp.SetNdivisions(5);
        h_jpt.SetNdivisions(5);

	h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; eta; N",80,-2,2)
	h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; phi; N",50,0,2*3.15)
	h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,200)
	h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,200)

#        print tin.GetEntriesFast()
	for i in range(tin.GetEntriesFast()):
            tin.GetEntry(i)
            tinmc.GetEntry(i)

            for n in range(len(tin.jpt)):
                if tin.jisleptag[n]==0 and abs(tin.jeta[n])<1.1:
                    for k in range(len(tinmc.jpt)):
                        dr = math.sqrt( (tin.jphi[n]-tinmc.jphi[k])*(tin.jphi[n]-tinmc.jphi[k]) + (tin.jeta[n]-tinmc.jeta[k])*(tin.jeta[n]-tinmc.jeta[k]) )
                        if dr<0.1:
                            h_jeratio.Fill( tin.je[n]/tinmc.je[k])
                    h_je.Fill( tin.je[n] )
                    h_jpt.Fill( tin.jpt[n] )
                    h_jp.Fill( tin.jp[n] )				
                    h_jeta.Fill( tin.jeta[n] )
                    h_jphi.Fill( tin.jphi[n] )
                    h_jmass.Fill( tin.jmass[n] )
                    h_jmass_sd.Fill( tin.jmass_sd[n] )	

#        print h_je.GetEntries()
        hists = []
        hists.append( h_jeratio )
	hists.append( h_je )
	hists.append( h_jpt )
#	hists.append( h_jp )
	hists.append( h_jeta )
	hists.append( h_jphi )
	hists.append( h_jmass )
	hists.append( h_jmass_sd )
	return hists


def MCinfo(tin):

    h_mzp = ROOT.TH1F("h_mh2","; H2 mass; N",100,0,10000);
    h_mzp.SetNdivisions(5);
    for i in range(tin.GetEntriesFast()):
        tin.GetEntry(i);
        h_mzp.Fill(tin.gen_mZp);
    return h_mzp;


## tin, PF, gen

def GetResolutions(tin,tin1,tin2,pf):

	h_jpres = ROOT.TH1F("h_jpres_"+pf, "; (P - PGEN)/PGEN; au", 50,-1,1)

	for i in range(tin.GetEntriesFast()):
		tin1.GetEntry(i)
		tin2.GetEntry(i)
                tin.GetEntry(i)
		for j in range(len(tin1.jpt)):
                    if tin1.jisleptag[j] == 0: 
                        for k in range(len(tin2.jpt)):
                            if tin2.jisleptag[k] == 0:
                                dr = math.sqrt( (tin1.jphi[j]-tin2.jphi[k])*(tin1.jphi[j]-tin2.jphi[k]) + (tin1.jeta[j]-tin2.jeta[k])*(tin1.jeta[j]-tin2.jeta[k]) )
                                if dr < 0.01: 
                                    for n in range(len(tin.jpt)):
                                        drpf = math.sqrt( (tin1.jphi[j]-tin.jphi[n])*(tin1.jphi[j]-tin.jphi[n]) + (tin1.jeta[j]-tin.jeta[n])*(tin1.jeta[j]-tin.jeta[n]) )
                                        if(drpf < 0.4 and tin.jisleptag[n]==0):
                                            h_jpres.Fill( (tin.jp[n] - tin2.jp[k])/tin2.jp[k] )				
                                                    
        return h_jpres

if __name__ == '__main__':

        directory = 'plots_radius' + jetRadius
        if not os.path.exists(directory):
                os.makedirs(directory)
        else:
                print 'the directory already exists! rememebr to clean up your work area'
                quit()

        fa = ROOT.TFile(inputFolder+"/radius"+jetRadius+"_of_PanPFA.root")
        tg = fa.Get("tGEN")
        tg_charged = fa.Get("tGEN_charged")
        tg_nonu    = fa.Get("tGEN_nonu")
        tg_response = fa.Get("tGEN_response")
        tg_resolution = fa.Get("tGEN_resolution")
        ta = fa.Get("tPFA")
        tc = fa.Get("tcalo")
        tt = fa.Get("ttrack")
        tmc = fa.Get("tMC")

        print "Getting hists..."
        hg = getHists(tg,tg,'_gen')
        hg_charged = getHists(tg_charged,tg,'_gen_charged')
        hg_nonu    = getHists(tg_nonu,tg,'_gen_nonu')
        hg_response = getHists(tg_response, tg,'_gen_response')
        hg_resolution = getHists(tg_resolution, tg,'_gen_resolution')

        ha = getHists(ta,tg,'_PF')
        ht = getHists(tt,tg,'_trk')	
        hc = getHists(tc,tg,'_calo')

        for i in range(len(ha)):
#            makeCanvas( [hg[i],hg_charged[i],hg_nonu[i],ha[i],ht[i], hc[i]], ['gen','gen_charged','gen_nonu','pf','track (associated) ','calo (associated)'] )
            makeCanvas( [hg_nonu[i],hg_response[i],hg_charged[i],hg_resolution[i], hc[i]], 

                        ['gen_nonu','gen_nonu_response','gen_charged (+ #gamma)','gen_charged_response (+ #gamma)','calo (associated)'] )

        hmc = MCinfo(tmc);
#        makeCanvas([hmc],['particlelevel'])

#        hres_pf    = GetResolutions(ta,ta,tg, "pf" )
#        hres_track = GetResolutions(tt,ta,tg, "track")
#        hres_cal   = GetResolutions(tc,ta,tg, "calo" )
#        makeCanvas( [hres_pf, hres_track, hres_cal], ['pf','track (associated)', 'calo (associated)'] )





