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
import sys

# make list of files..
import glob

for arg in sys.argv: 
    print arg



filename = sys.argv[1]
print filename
outputtemp=filename
directory=outputtemp.rsplit('/', 1)[1]
print directory

files = glob.glob(filename+"/*.slcio")

factory = LCFactory.getInstance()
reader = factory.createLCReader()
if not os.path.exists(directory):
    os.makedirs(directory)
else:
    print 'the directory already exists! rememebr to clean up your work area'
    quit()

fileOutStatus1         = open(directory+'/of_status1.dat','w');
fileOutStatus3         = open(directory+'/of_status3.dat','w');
fileOutPanPFA          = open(directory+'/of_PanPFA.dat','w');
fileOutPanPFA_Tracks   = open(directory+'/of_PanPFA_Tracks.dat','w');
fileOutPanPFA_Calo     = open(directory+'/of_PanPFA_Calo.dat','w');
fileOutClusters        = open(directory+'/of_CaloClusters.dat','w');
fileOutTracks          = open(directory+'/of_Tracks.dat','w');
fileOutECaloR          = open(directory+'/of_ECaloHits_r.dat','w');
fileOutHCaloR          = open(directory+'/of_HCaloHits_r.dat','w');



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
		col    = evt.getCollection("MCParticle")
		colPF  = evt.getCollection("PandoraPFOCollection");
		colCl  = evt.getCollection("ReconClusters");
		colTr  = evt.getCollection("Tracks");
		colECB = evt.getCollection("EcalBarrelHits");
                colHCB = evt.getCollection("HcalBarrelHits");

		nMc=col.getNumberOfElements();
		nPF=colPF.getNumberOfElements();
		nCl=colCl.getNumberOfElements();
		nTr=colTr.getNumberOfElements();                
		nECB=colECB.getNumberOfElements();		
		nHCB=colHCB.getNumberOfElements();

#		print " ----------- ------------ ",nEvent
#		print "number of particles = ", nMc, nPF, nTr
		fileOutStatus1.write( '-99 0 0 0 0 \n' );
		fileOutStatus3.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA_Calo.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA_Tracks.write( '-99 0 0 0 0 \n' );
		fileOutClusters.write( '-99 0 0 0 0 \n' );
		fileOutTracks.write( '-99 0 0 0 0 \n' );
                fileOutECaloR.write( '-99 0 0 0 0 \n' ); 
                fileOutHCaloR.write( '-99 0 0 0 0 \n' ); 

		for i in range(nMc): # loop over all particles 
			par=col.getElementAt(i)
			if par.getGeneratorStatus() == 1: 
                            fileOutStatus1.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(par.getPDG(),par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy()) )
		
		for i in range(nPF): # loop over all particles
			parPF=colPF.getElementAt(i) 
			fileOutPanPFA.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getType(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy()) )

			## Look for CaloClusters associated with PF			
			if parPF.getClusters().size() > 1:
                                print "clusters: ", nEvent, i, parPF.getClusters().size()                                
			if parPF.getClusters().size() > 0: 
                            fileOutPanPFA_Calo.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getClusters()[0].getType(),parPF.getClusters()[0].getPosition()[0],parPF.getClusters()[0].getPosition()[1],parPF.getClusters()[0].getPosition()[2],parPF.getClusters()[0].getEnergy()) )

                        ## Look for tracks associated with PF 
			if parPF.getTracks().size() > 1:
                                print nEvent, i, parPF.getTracks().size()                                
			if parPF.getTracks().size() > 0:
                                trackpt  = math.fabs( (3e-4*5/parPF.getTracks()[0].getOmega()));
                                cotTheta = parPF.getTracks()[0].getTanLambda();
                                trackp   = trackpt * math.sqrt(1 + cotTheta*cotTheta );
                                phi      = parPF.getTracks()[0].getPhi();
                                trackpx  = trackpt*math.cos(phi);
                                trackpy  = trackpt*math.sin(phi);
                                trackpz  = trackpt*cotTheta;
                                trackE   = math.sqrt(trackp*trackp + 0.13957*0.13957)

                                fileOutPanPFA_Tracks.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getTracks()[0].getType(),trackpx,trackpy,trackpz,trackE))

                # print out calorimeter hits
		for i in range(nCl):
			parCl = colCl.getElementAt(i);
			fileOutClusters.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parCl.getType(),parCl.getPosition()[0],parCl.getPosition()[1],parCl.getPosition()[2],parCl.getEnergy()) )

#		for i in range(nECB):
#			parCl = colECB.getElementAt(i);
#			fileOutECaloR.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(0,parCl.getPosition()[0],parCl.getPosition()[1],parCl.getPosition()[2],parCl.getEnergy()) )

#		for i in range(nHCB):
#			parCl = colHCB.getElementAt(i);
#			fileOutHCaloR.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(1,parCl.getPosition()[0],parCl.getPosition()[1],parCl.getPosition()[2],parCl.getEnergy()) )


#		if nEvent > 10: break;
reader.close() # close the file



