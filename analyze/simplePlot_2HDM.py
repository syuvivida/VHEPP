import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15);
ROOT.gStyle.SetPalette(1);
import math

def makeCanvas(hists, tags):

	colors = [1,2,4,5,6,7]

	leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	tmax = -999;
	for i in range(len(hists)):
		leg.AddEntry(hists[i],tags[i],'l');
		if hists[i].GetMaximum() > tmax: tmax = hists[i].GetMaximum();

	c = ROOT.TCanvas("c","c",1000,800);
	hists[0].SetMinimum(0);
	hists[0].SetMaximum(tmax*1.2);
	for i in range(len(hists)):
		if i == 0: hists[i].Draw();
		else: 
			hists[i].SetLineColor( colors[i] );
			hists[i].Draw('sames');

	leg.Draw();
	c.SaveAs("plots/"+hists[0].GetName()+".pdf")
	c.SaveAs("plots/"+hists[0].GetName()+".png")

def getHists(tin,postfix):

	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,6000);
	h_jpt.SetNdivisions(5);
	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,0,6000);
	h_jp.SetNdivisions(5);
	h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; eta; N",40,-1,1);
	h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; phi; N",50,0,2*3.15);
	h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,200);
#	h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,200);

	# print tin.GetEntriesFast();
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i);
		for j in range(len(tin.jpt)):
			# print tin.jisleptag[j]
			if tin.jisleptag[j] == 0: 
				# print tin.jpt[j]
				h_jpt.Fill( tin.jpt[j] );
				h_jp.Fill( tin.jp[j] );				
				h_jeta.Fill( tin.jeta[j] );
				h_jphi.Fill( tin.jphi[j] );
				h_jmass.Fill( tin.jmass[j] );
#				h_jmass_sd.Fill( tin.jmass_sd[j] );	

	hists = [];
	hists.append( h_jpt );
	hists.append( h_jp );
	hists.append( h_jeta );
	hists.append( h_jphi );
	hists.append( h_jmass );
#	hists.append( h_jmass_sd );
	return hists;

def MCinfo(tin):

	h_mzp = ROOT.TH1F("h_mh2","; H2 mass; N",100,4000,6000);
	h_mzp.SetNdivisions(5);
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i);
		h_mzp.Fill(tin.gen_mZp);
	return h_mzp;

def GetResolutions(tin1,tin2,pf):

	h_jpres = ROOT.TH1F("h_jpres_"+pf, "; (P - PGEN)/PGEN; au", 50,-0.5,0.5);

	for i in range(tin1.GetEntriesFast()):
		tin1.GetEntry(i);
		tin2.GetEntry(i);

		for j in range(len(tin1.jpt)):
			if tin1.jisleptag[j] == 0: 
				for k in range(len(tin2.jpt)):

					dr = math.sqrt( (tin1.jphi[j]-tin2.jphi[k])*(tin1.jphi[j]-tin2.jphi[k]) + (tin1.jeta[j]-tin2.jeta[k])*(tin1.jeta[j]-tin2.jeta[k]) );
					if dr < 0.4: 
						h_jpres.Fill( (tin1.jp[j] - tin2.jp[k])/tin2.jp[k] );				
				
	return h_jpres;

if __name__ == '__main__':

	fa = ROOT.TFile("dat/of_PanPFA.root");
	ta = fa.Get("tPFA");
	tb = fa.Get("tGEN");
	tc = fa.Get("tcalo");
	tmc = fa.Get("tMC");


	print "Getting hists..."
	ha = getHists(ta,'a');
	hb = getHists(tb,'b');
	hc = getHists(tc,'c');	

	for i in range(len(ha)):
		makeCanvas( [ha[i],hb[i],hc[i]], ['pf','gen','calo'] )

	hmc = MCinfo(tmc);
	makeCanvas([hmc],['particlelevel'])
	# hres = Resolutions(tb,ta,tc);

	hres_pf = GetResolutions( ta, tb, "pf" );
	hres_cal = GetResolutions( tc, tb, "calo" );
	makeCanvas( [hres_pf,hres_cal], ['pf','calo'] )




