import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15);
ROOT.gStyle.SetPalette(1);

def readEvent(fin):
	
	event = [];
	for line in fin:
		if '####' in line: break;
		else:
			linearray = line.strip().split();
			event.append( linearray );
	
	return event

def makeEvent(listofparticles,tag):

	ev = ROOT.TH2F("ev",";#eta;#phi",200,-5,5,126,-3.15,+3.15);

	for p in listofparticles:
		fv = ROOT.TLorentzVector();
		fv.SetPxPyPzE(float(p[1]),float(p[2]),float(p[3]),float(p[4]));
		# print "eta,phi = ", fv.Eta(), fv.Phi(), fv.Energy();
		# print ev.FindBin(fv.Eta(),fv.Phi());
		curbincontent = ev.GetBinContent(ev.FindBin(fv.Eta(),fv.Phi()));
		ev.SetBinContent( ev.FindBin(fv.Eta(),fv.Phi()), fv.Energy()+curbincontent );

	# for i in range(ev.GetNbinsX()):
	# 	for j in range(ev.GetNbinsY()):
	# 		print i,j,ev.GetBinContent(i+1,j+1);

	# print ev.GetMaximum()

	can = ROOT.TCanvas('can','can',1200,800);
	# ev.SetMinimum(-1.);
	# ev.SetMarkerStyle(21);
	ev.SetMaximum(6000.);
	ev.Draw('colz');
	ROOT.gPad.SetLogz();
	can.SaveAs("displays/display_"+tag+".pdf");

	del ev;
	del can;

if __name__ == '__main__':

	if_mc3 = open('of_status3.dat','r');
	# if_mc3 = open('of_PanPFA.dat','r');
	tag = 'mc3';
	# eventDisplay

	eventnumber = 0;
	while(1):

		event_mc3 = readEvent(if_mc3);
		# print event_mc3;
		if len(event_mc3) == 0 and eventnumber != 0: break;

		print "--------------------";
		print "event ",eventnumber,len(event_mc3);
		# for p in event_mc3:
		# 	print p
		
		if eventnumber > 0: makeEvent(event_mc3,str(eventnumber)+"_"+tag);

		eventnumber += 1;





	test = ROOT.TH2F("test","test",2,0,2,2,0,2);
	test.SetBinContent(1,1,5)
	test.SetBinContent(1,2,25)
	test.SetBinContent(2,1,10)
	test.SetBinContent(2,2,15)
	cant = ROOT.TCanvas("cant","cant",800,800);
	test.Draw('colz');
	cant.SaveAs("cant.pdf");


