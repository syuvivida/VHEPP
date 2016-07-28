
import ROOT
import tdrstyle
import math
import os

import sys
for arg in sys.argv: 
    print arg

inputFolder = sys.argv[1]
jetRadius = sys.argv[2]
sample    = sys.argv[3]

tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gStyle.SetPalette(1)

## Obtain proper dr
def myDeltaR(eta1, eta2, phi1, phi2):

   deta = eta1-eta2;
   dphi = phi1-phi2;

   if dphi >= math.pi: 
       dphi = dphi-math.pi*2;
   if dphi < -math.pi:
       dphi = dphi+math.pi*2;

   dr = math.sqrt(deta*deta+dphi*dphi)
   return dr


## Obtain FWHM

def FWHM(hist):
    bin1 = hist.FindFirstBinAbove(hist.GetMaximum()/2);
    bin2 = hist.FindLastBinAbove(hist.GetMaximum()/2);
    fwhm = hist.GetBinCenter(bin2) - hist.GetBinCenter(bin1);
    return fwhm

## Making final figures and setup final display

def makeCanvas(sample, hists, tags):
    colors = [1,2,4,8,6,7,8]

#    if hists[0].GetName() == "h_jp_gen_nonu" or hists[0].GetName() == "h_jpt_gen_nonu" or hists[0].GetName() == "h_je_gen_nonu":            
#        leg = ROOT.TLegend(0.4,0.6,0.8,0.9)
#    else:
    if hists[0].GetName() == "h_jeratio_gen_response":
        leg = ROOT.TLegend(0.18,0.4,0.8,0.9)
    elif ((sample =='1' or sample =='2' or sample =='3' or sample =='4') and (hists[0].GetName() == "h_jp_gen_nonu" or hists[0].GetName() == "h_jpt_gen_nonu" or hists[0].GetName() == "h_je_gen_nonu")) or hists[0].GetName() == "h_jmass_gen_nonu" or hists[0].GetName() == "h_jmass_sd_gen_nonu":
        leg = ROOT.TLegend(0.45,0.6,0.85,0.9)
    else:
        leg = ROOT.TLegend(0.18,0.6,0.58,0.9)
    leg.SetHeader("Jet radius="+jetRadius)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    tmax = -999


    for i in range(len(hists)):
        if hists[0].GetName() == "h_jeratio_gen_response":
            tagname="RMS = " + "{:.3f}".format(hists[i].GetRMS())
            tagname2="FWHM = " + "{:.3f}".format(FWHM(hists[i]))
            leg.AddEntry(hists[i],tags[i],'l')
            leg.AddEntry(None, tagname,'')
            leg.AddEntry(None, tagname2,'')
        elif hists[0].GetName() != "h_jeratio_gen_nonu" :
            leg.AddEntry(hists[i],tags[i],'l')
        if hists[i].GetMaximum() > tmax:
            tmax = hists[i].GetMaximum()
            
    if hists[0].GetName() == "h_jeta_gen_nonu" or hists[0].GetName() == "h_jphi_gen_nonu": 
        hists[0].SetMaximum(tmax*2.0)
#    elif hists[0].GetName() == "h_jeratio_gen_nonu":
#        hists[1].SetMaximum(tmax*1.2)
    else:
        hists[0].SetMaximum(tmax*1.2)
            
    c = ROOT.TCanvas("c","c",1000,800)
    for i in range(len(hists)):
        hists[i].SetMinimum(0)
        hists[i].GetYaxis().SetDecimals()
        hists[i].SetLineWidth(3)
        hists[i].SetLineStyle(1)
        if i == 0 or i==2:
            hists[i].SetLineStyle(2+i)
        if i == 0:
            hists[i].DrawNormalized('hist')
        elif hists[0].GetName() == "h_jeratio_gen_nonu" and i ==1:
            hists[i].DrawNormalized('hist')            
        else: 
            hists[i].SetLineColor( colors[i] )
            hists[i].DrawNormalized('histsames')

    leg.Draw()

    c.SaveAs(directory+"/radius"+jetRadius+"_"+hists[0].GetName()+".pdf")
    c.SaveAs(directory+"/radius"+jetRadius+"_"+hists[0].GetName()+".png")

## tin, PF, gen
### Filling histograms #########
def getHists(sample,tin,tinmc,postfix):


    titleratio = "; E_{jet}/E_{jet}^{true}; N"
    titleje    = "; jet energy [GeV]; N"
    titlejpt   = "; jet p_{T} [GeV]; N"
    titlejp    = "; jet momentum [GeV]; N"
    titlejmass = "; jet mass [GeV]; N"
    titlejsd   = "; soft drop mass (#beta = 0) [GeV]; N"

    h_jeratio  = ROOT.TH1F("h_jeratio"+postfix, titleratio, 50,0.6,1.1)

# 2HDM
    if sample =='0':
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,1000,3000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,3000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,1000,3000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,180)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,180)
        
    elif sample =='6':
## Z' = 10 TeV sqrt(s)=10 TeV
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,2000,6000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,6000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,2000,6000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,800)
# Z' = 5 TeV sqrt(s)=40 TeV
    elif sample =='1':
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,0,10000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,10000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,0,10000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,800)
# Z' = 10 TeV sqrt(s)=40 TeV
    elif sample =='2':
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,0,20000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,20000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,0,20000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,800)
# Z' = 20 TeV sqrt(s)=40 TeV
    elif sample =='3':
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,0,35000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,35000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,0,35000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,800)
# Z' = 30 TeV sqrt(s)=40 TeV
    elif sample =='4':
        h_je       = ROOT.TH1F("h_je"+postfix,      titleje,    50,0,45000);
        h_jpt      = ROOT.TH1F("h_jpt"+postfix,     titlejpt,   50,0,45000)
        h_jp       = ROOT.TH1F("h_jp"+postfix,      titlejp,    50,0,45000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,   titlejmass, 50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,titlejsd,   50,0,800)
    else:
	h_je       = ROOT.TH1F("h_je"+postfix,"; jet energy; N",50,10000,22000);
	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",50,0,22000)
	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",50,10000,22000)
        h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,800)
        h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,800)

    
    h_jeratio.SetNdivisions(5);
    h_je.SetNdivisions(5);
    h_jp.SetNdivisions(5);
    h_jpt.SetNdivisions(5);
    
    h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; jet #eta; N",80,-2,2)
    h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; jet #phi; N",50,0,2*3.15)

#        print tin.GetEntriesFast()
    for i in range(tin.GetEntriesFast()):
        tin.GetEntry(i)
        tinmc.GetEntry(i)

        for n in range(len(tin.jpt)):
            if tin.jisleptag[n]==0 and abs(tin.jeta[n])<1.1:
                for k in range(len(tinmc.jpt)):
                    
                    dr = myDeltaR(tin.jeta[n],tinmc.jeta[k],tin.jphi[n],tinmc.jphi[k])
                    
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
### Resolution function ############
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
                        dr = myDeltaR(tin1.jeta[j],tin2.jeta[k],tin1.jphi[j],tin2.jphi[k])
                        if dr < 0.01: 
                            for n in range(len(tin.jpt)):
                                drpf = myDeltaR(tin1.jeta[j],tin.jeta[n],tin1.jphi[j],tin.jphi[n])
                                if(drpf < 0.4 and tin.jisleptag[n]==0):
                                    h_jpres.Fill( (tin.jp[n] - tin2.jp[k])/tin2.jp[k] )				
                                                    
    return h_jpres


########## Main Function #############
if __name__ == '__main__':

    directory = 'plots_radius' + jetRadius
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        print 'the directory already exists! rememebr to clean up your work area'
        quit()



    fa = ROOT.TFile(inputFolder+"/radius"+jetRadius+"_response.root")
    tg = fa.Get("tGEN")
    tg_nonu    = fa.Get("tGEN_nonu")
    tg_response = fa.Get("tGEN_response")
    tg_charged = fa.Get("tGEN_charged")
    tg_resolution = fa.Get("tGEN_resolution")
    ta = fa.Get("tPFA")
    tc = fa.Get("tcalo")
    tt = fa.Get("ttrack")
    tmc = fa.Get("tMC")



    ta = fa.Get("tPFA")
    tc = fa.Get("tcalo")
#    traw = fa.Get("trawhits")
#    traw2 = fa.Get("trawhits2")
    tt = fa.Get("ttrack")
    tmc = fa.Get("tMC")

    print "Getting hists..."
    hg = getHists(sample,tg,tg_nonu,'_gen')
    hg_nonu    = getHists(sample,tg_nonu,tg_nonu,'_gen_nonu')
    hg_response = getHists(sample,tg_response, tg_nonu,'_gen_response')
    hg_charged = getHists(sample,tg_charged,tg_nonu,'_gen_charged')
    hg_resolution = getHists(sample,tg_resolution, tg_nonu,'_gen_resolution')

    ha = getHists(sample,ta,tg_nonu,'_PF')
    ht = getHists(sample,tt,tg_nonu,'_trk')	
    hc = getHists(sample,tc,tg_nonu,'_calo')
#    hraw = getHists(sample,traw,tg_nonu,'_raw')
#    hraw2 = getHists(sample,traw2,tg_nonu,'_raw2')

    for i in range(len(ha)):
#        makeCanvas( sample, [hg_nonu[i],hg_response[i],hc[i], hraw[i], hraw2[i]],                         
 #                   ['gen (no #nu)','gen response (no #nu)','calo clusters (PF)', 'calo hits (fastjet)', 'calo hits (simple)' ] )
#        makeCanvas( sample, [hg_nonu[i],hg_response[i],hc[i], hraw2[i]],                         
 #                   ['gen (no #nu)','gen response (no #nu)','calo clusters (PF)', 'calo hits (simple)' ] )
        makeCanvas( sample, [hg_response[i],hg_charged[i],hg_resolution[i], hc[i]],                         
                      ['gen response (no #nu)','gen (charged + #gamma)','gen response (charged+ #gamma)','PF calo'] )

    hmc = MCinfo(tmc);
#        makeCanvas([hmc],['particlelevel'])




