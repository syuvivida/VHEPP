# Example: read truth-level particles and create Z peak  
# Use http://atlaswww.hep.anl.gov/hepsim/show.php?item=146&reco=rfull001 
#  Get 4 files in 2 threads (see HepSim manual):
#  > hs-get gev250ee_pythia6_zpole_ee%rfull001 gev250ee_pythia6_zpole_ee 2 4 
# S.Chekanov (ANL)

from org.lcsim.lcio import LCIOReader
from hep.io.sio import SIOReader
from hep.lcio.implementation.sio import SIOLCReader
from hep.lcio.implementation.io import LCFactory
from hep.lcio.event import * 
from hep.lcio.io import *
from jhplot import *  # import graphics
from hephysics.particle import LParticle
import math
import os

# make list of files..
import glob
files=glob.glob("/home/syu/muon_collider/single_pion_analog/rfull005/pgun_pi1000gev*.slcio")
factory = LCFactory.getInstance()
reader = factory.createLCReader()
directory = 'dat'
if not os.path.exists(directory):
    os.makedirs(directory)
else:
	print 'the directory already exists! rememebr to clean up your work area'
	quit()

fileOutStatus1  = open('dat/of_status1.dat','w');
fileOutStatus3  = open('dat/of_status3.dat','w');
fileOutPanPFA   = open('dat/of_PanPFA.dat','w');
fileOutPanPFA_Tracks   = open('dat/of_PanPFA_Tracks.dat','w');
fileOutPanPFA_CaloHits   = open('dat/of_PanPFA_CaloHits.dat','w');
fileOutClusters = open('dat/of_CaloHits.dat','w');
fileOutTracks   = open('dat/of_Tracks.dat','w');


nEvent=0
for f in files:
	print "Open file=",f
	reader.open(f) 
	while(1):
		evt=reader.readNextEvent()
		if (evt == None): break
		nEvent=nEvent+1
		# print " file event: ",evt.getEventNumber(), " run=",evt.getRunNumber()
#		if (nEvent%100==0): print "# Event: ",nEvent
		col = evt.getCollection("MCParticle")
		colPF = evt.getCollection("PandoraPFOCollection");
		colCl = evt.getCollection("ReconClusters");
		colTr = evt.getCollection("Tracks");
		nMc=col.getNumberOfElements()
		nPF=colPF.getNumberOfElements();
		nCl=colCl.getNumberOfElements();
		nTr=colTr.getNumberOfElements();
#		print " ----------- ------------ ",nEvent
#		print "number of particles = ", nMc, nPF, nTr
		fileOutStatus1.write( '-99 0 0 0 0 \n' );
		fileOutStatus3.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA_CaloHits.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA_Tracks.write( '-99 0 0 0 0 \n' );
		fileOutClusters.write( '-99 0 0 0 0 \n' );
		fileOutTracks.write( '-99 0 0 0 0 \n' );
		for i in range(nMc): # loop over all particles 
			par=col.getElementAt(i)
			if par.getGeneratorStatus() == 3: 
				pdg=par.getPDG();
				charge=par.getCharge();
				if par.getParents().size() > 0:
					if math.fabs(par.getParents()[0].getPDG()) == 24 or math.fabs(par.getPDG()) == 24 or par.getPDG() > -999999:
						print i, "pdg=",pdg, ", status=",par.getGeneratorStatus(), " mass=",par.getMass(), " mother=",par.getParents()[0].getPDG(), " momsize=",par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy(); 
						fileOutStatus3.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(par.getPDG(),par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy()) )
			if par.getGeneratorStatus() == 1: 
				if (math.fabs(par.getPDG()) == 12 or math.fabs(par.getPDG()) == 14 or math.fabs(par.getPDG()) == 16): continue;
				fileOutStatus1.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(par.getPDG(),par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy()) )
		
		for i in range(nPF): # loop over all particles
			parPF=colPF.getElementAt(i) 
			fileOutPanPFA.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getType(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy()) )
			## Look for CaloClusters associated with PF			
			ECluster = -99;
			if parPF.getClusters().size() > 0: 
				for j in range(parPF.getClusters()[0].getCalorimeterHits().size()):
					fileOutPanPFA_CaloHits.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getClusters()[0].getCalorimeterHits()[j].getType(),parPF.getClusters()[0].getCalorimeterHits()[j].getPosition()[0],parPF.getClusters()[0].getCalorimeterHits()[j].getPosition()[1],parPF.getClusters()[0].getCalorimeterHits()[j].getPosition()[2],parPF.getClusters()[0].getCalorimeterHits()[j].getEnergy()) )
			if parPF.getTracks().size() > 1:
				print i, parPF.getTracks().size()

#			if parPF.getTracks().size() > 0:
#				for k in range(parPF.getTracks().size()):
#					fileOutPanPFA_Tracks.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getTracks()[k].getType(),parPF.getTracks()[k].getMomentum()[0],parPF.getTracks()[k].getMomentum()[1],parPF.getTracks()[k].getMomentum()[2],parPF.getTracks()[k].getEnergy()))
					
		for i in range(nCl):
			parCl = colCl.getElementAt(i);
			# print "position = ",len(parCl.getPosition()), parCl.getPosition()[0], parCl.getPosition()[1], parCl.getPosition()[2], parCl.getEnergy();
			# perp2 = parCl.getPosition()[0]*parCl.getPosition()[0] + parCl.getPosition()[1]*parCl.getPosition()[1]
			# print "phi = ", math.atan2(parCl.getPosition()[1],parCl.getPosition()[0]), "and theta = ", math.atan2(math.sqrt(perp2),parCl.getPosition()[2])
			fileOutClusters.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parCl.getType(),parCl.getPosition()[0],parCl.getPosition()[1],parCl.getPosition()[2],parCl.getEnergy()) )

		# for i in range(nTr):
		# 	parTr = colTr.getElementAt(i);
		# 	fileOutTracks.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parTr.getType(),parTr.getMomentum()[0],parTr.getMomentum()[1],parTr.getMomentum()[2],parTr.getEnergy()) )

#		if nEvent > 1000: break;
reader.close() # close the file



