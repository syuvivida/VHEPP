import ROOT
import tdrstyle
import math
import os

tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gStyle.SetPalette(1)

def makeCanvas(hists, tags):
    colors = [1,2,4,5,6,7]
            
    leg = ROOT.TLegend(0.2,0.7,0.6,0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    tmax = -999

    for i in range(len(hists)):
        leg.AddEntry(hists[i],tags[i],'l')
        if hists[i].GetMaximum() > tmax:
            tmax = hists[i].GetMaximum()
            
    c = ROOT.TCanvas("c","c",1000,800)
    for i in range(len(hists)):
        if i == 0:
            hists[i].Draw()
        else: 
            hists[i].SetLineColor( colors[i] )
            hists[i].Draw('sames')

        hists[0].Draw('sames')
        leg.Draw()

        hists[0].SetMinimum(0)

        if hists[0].GetName() == "h_jeta_gen" or hists[0].GetName() == "h_jpres_pf" :
            hists[0].SetMaximum(tmax*1.5)
        else:
            hists[0].SetMaximum(tmax*1.3)

        c.SaveAs("plots/"+hists[0].GetName()+".pdf")
        c.SaveAs("plots/"+hists[0].GetName()+".png")

def getHists(tin,tin1,tin2,postfix):
### For 1 GeV single pion
#	h_je      = ROOT.TH1F("h_je"+postfix,"; jet energy; N",100,0,1.5)
#	h_je.SetNdivisions(5)
#	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,1.5)
#	h_jpt.SetNdivisions(5)
#	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,0,1.5)
#	h_jp.SetNdivisions(5)


### For 10 GeV single pion
	h_je      = ROOT.TH1F("h_je"+postfix,"; jet energy; N",150,0,15);
	h_je.SetNdivisions(5)
	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",150,0,15)
	h_jpt.SetNdivisions(5)
	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",150,0,15)
	h_jp.SetNdivisions(5)

        

	h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; eta; N",80,-2,2)
	h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; phi; N",50,0,2*3.15)
	h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,20)
	h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,200)

#        print tin.GetEntriesFast()
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i)
                tin1.GetEntry(i)
                tin2.GetEntry(i)
#                print len(tin1.jpt), len(tin2.jpt)
                for j in range(len(tin1.jpt)):
#                    print "j=", j
                    for k in range(len(tin2.jpt)):
#                        print "k=", k
                        dr = math.sqrt( (tin1.jphi[j]-tin2.jphi[k])*(tin1.jphi[j]-tin2.jphi[k]) + (tin1.jeta[j]-tin2.jeta[k])*(tin1.jeta[j]-tin2.jeta[k]) )
#                        print "dr=", dr
                        if dr < 0.01: 
                            for n in range(len(tin.jpt)):
                                h_je.Fill( tin.je[n] )
                                h_jpt.Fill( tin.jpt[n] )
                                h_jp.Fill( tin.jp[n] )				
                                h_jeta.Fill( tin.jeta[n] )
                                h_jphi.Fill( tin.jphi[n] )
                                h_jmass.Fill( tin.jmass[n] )
                                h_jmass_sd.Fill( tin.jmass_sd[n] )	

        print h_je.GetEntries()
        hists = []
	hists.append( h_je )
	hists.append( h_jpt )
	hists.append( h_jp )
	hists.append( h_jeta )
	hists.append( h_jphi )
	hists.append( h_jmass )
	hists.append( h_jmass_sd )
	return hists


def GetResolutions(tin1,tin2,pf):

	h_jpres = ROOT.TH1F("h_jpres_"+pf, "; (P - PGEN)/PGEN; au", 50,-0.2,0.2)

	for i in range(tin1.GetEntriesFast()):
		tin1.GetEntry(i)
		tin2.GetEntry(i)

		for j in range(len(tin1.jpt)):
			if tin1.jisleptag[j] == 0: 
				for k in range(len(tin2.jpt)):

					dr = math.sqrt( (tin1.jphi[j]-tin2.jphi[k])*(tin1.jphi[j]-tin2.jphi[k]) + (tin1.jeta[j]-tin2.jeta[k])*(tin1.jeta[j]-tin2.jeta[k]) )
					if dr < 0.01: 
						h_jpres.Fill( (tin1.jp[j] - tin2.jp[k])/tin2.jp[k] )				
				
	return h_jpres

if __name__ == '__main__':

        directory = 'plots'
        if not os.path.exists(directory):
                os.makedirs(directory)
        else:
                print 'the directory already exists! rememebr to clean up your work area'
                quit()

        fa = ROOT.TFile("dat/of_PanPFA.root")
        ta = fa.Get("tPFA")
        tg = fa.Get("tGEN")
        tc = fa.Get("tcalo")
        tt = fa.Get("ttrack")
        tmc = fa.Get("tMC")

        print "Getting hists..."
        ha = getHists(ta,ta,tg,'_PF')
        hg = getHists(tg,ta,tg,'_gen')
        hc = getHists(tc,ta,tg,'_calo')
        ht = getHists(tt,ta,tg,'_trk')	

        for i in range(len(ha)):
            makeCanvas( [hg[i],ha[i],ht[i], hc[i]], ['gen','pf','track (associated) ','calo (associated)'] )


        hres_pf = GetResolutions( ta, tg, "pf" )
        hres_track = GetResolutions( tt, tg, "track")
        hres_cal = GetResolutions( tc, tg, "calo" )
        makeCanvas( [hres_pf, hres_track, hres_cal], ['pf','track (associated)', 'calo (associated)'] )




