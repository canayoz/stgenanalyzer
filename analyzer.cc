//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Wed Jul  9 15:44:25 2014 by mkanalyzer.py
// Author:      Halil Gamsizkan
//-----------------------------------------------------------------------------
#include "analyzer.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line

  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples

  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

  // Select variables to be read

  selectVariables(stream);


  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the
  // following line

  TApplication app("analyzer", &argc, argv);

  /**
	 Notes 1
	 -------
	 1. Use
	   ofile = outputFile(cmdline.outputfile, stream)

	 to skim events to output file in addition to writing out histograms.

	 2. Use
	   ofile.addEvent(event-weight)

	 to specify that the current event is to be added to the output file.
	 If omitted, the event-weight is defaulted to 1.

	 3. Use
		ofile.count(cut-name, event-weight)

	 to keep track, in the count histogram, of the number of events
	 passing a given cut. If omitted, the event-weight is taken to be 1.
	 If you want the counts in the count histogram to appear in a given
	 order, specify the order, before entering the event loop, as in
	 the example below

		ofile.count("NoCuts", 0)
		ofile.count("GoodEvent", 0)
		ofile.count("Vertex", 0)
		ofile.count("MET", 0)

     Notes 2
	 -------
	 By default all variables are saved. Before the event loop, you can use
  
       select(objectname)
	  
     e.g.,
	
       select("GenParticle")
  
     to declare that you intend to select objects of this type. The
	 selection is done using

       select(objectname, index)
	  
     e.g.,
	  
       select("GenParticle", 3),
  
     which is called within the event loop. Call saveSelectedObjects()
	 before a call to addEvent if you wish to save the selected objects.
	 All other objects are saved by default.
	 
	 NB: If you declare your intention to select objects of a given type
	     by calling select(objectname), but subsequently fail to select
	     them using select(objectname, index) then none will be saved!
  */

  outputFile ofile(cmdline.outputfilename);

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------


using namespace std;

string pdg[25000];
for (int q=0; q<25000; q++){pdg[q]="?";}
pdg[2]="u";
pdg[1]="d";
pdg[3]="s";
pdg[4]="c";
pdg[5]="b";
pdg[6]="t";
pdg[11]="e-";
pdg[12]="nue";
pdg[13]="mu-";
pdg[14]="numu";
pdg[15]="tau-";
pdg[16]="nutau";
pdg[21]="g";
pdg[22]="gamma";
pdg[23]="Z";
pdg[24]="W+";
  //nevents=2;
  for(int entry=0; entry < nevents; ++entry)
	{
	  // Read event into memory
	  stream.read(entry);

	  // NB: call to clear object selection map (indexmap)
	  initialize();
	  
	  // Uncomment the following line if you wish to copy variables into
	  // structs. See the header file analyzer.h to find out what structs
	  // are available. Alternatively, you can call individual fill functions.
	 fillObjects();
	 
	 
	 cout << "Event #" << entry << ", #GJ:" << GenJet.size() << ", #GenPart.s:" << GenParticle.size() << endl;
	 //cout << "New event\n";                                                                                                                     
          // for (int q=0; q<GenJet.size(); q++) {                                                                                                   
                    //  std::cout << GenJet[q] << std::endl; 
          int myT=0;
	  int myj=0;
	  int mybbar=0;
	  int lastT=0;
	  int mytb=0;
	  int mytW=0;
	  int myWl=0;
	  int myWnu=0;
	  for (int q=0; q<GenParticle.size() ; q++) {
	    if (GenParticle[q].pdgId==21 || GenParticle[q].pdgId==22) continue; // gluons & photons
	    if (q<10 && GenParticle[q].numberOfMothers == 2 && GenParticle[q].pdgId == 6) myT=q; // found initial top
	    if (q<10 && GenParticle[q].numberOfMothers == 2 && GenParticle[q].pdgId == -5) mybbar=q; // found my bbar
	    if (q<10 && GenParticle[q].numberOfMothers == 2 && abs(GenParticle[q].pdgId) < 5) myj=q; // found my tagging jet
	    
	    if (  myT>0
		  && (GenParticle[q].pdgId == 5 && mytb == 0)
		  && (mytW == 0 && GenParticle[q+1].pdgId==24)
	       ) {
	      mytb=q; // found the top decay products my b (t->Wb)
	      mytW=q+1; // W is always next
	    }
	    //if (myT>0 && GenParticle[q].pdgId == 24 && mytW == 0) mytW=q; // found my W+ (t->Wb)
	    
	    // 
	    if (  mytW > 0
		  && ((GenParticle[q].pdgId == -11 || GenParticle[q].pdgId == -13  || GenParticle[q].pdgId == -15 ) && myWl == 0)
		  && ((GenParticle[q+1].pdgId == 12 || GenParticle[q+1].pdgId == 14 || GenParticle[q+1].pdgId == 16 ) && myWnu == 0)
	       ) {
	      myWl=q; // found my l+ (W->l+nu)
	      myWnu=q+1; // neutrino is always next
	    }
	    //if (mytW > 0 && (GenParticle[q].pdgId == 12 || GenParticle[q].pdgId == 14 ) && myWnu == 0) myWnu=q; // found my nu (W->l+nu)
	    //
	    if (GenParticle[q].pdgId == 6) lastT=q; // mark the last occurence of t quark
	    
/*            if (GenParticle[q].pdgId == 6) lastT=q; // mark the last occurence of t quark
			std::cout << q << "   pdg:" << GenParticle[q].pdgId << " (" << pdg[abs(GenParticle[q].pdgId)] << ")" 
                              << " - pt:" << GenParticle[q].pt
                              << " - eta:" << GenParticle[q].eta
			      << " - #mothers:" << GenParticle[q].numberOfMothers
			      << " - top:" << myT << " bbar:" << mybbar << " tjet:" << myj
			      << " t-b:" << mytb << " tW:" << mytW << " Wl:" << myWl << " Wnu" << myWnu
                             << std::endl;
*/
          }

	  cout << "Decay information:\n"
	       << "Top..........: " << myT << "   pdg:" << GenParticle[myT].pdgId << " (" << pdg[abs(GenParticle[myT].pdgId)] << ")" << " - pt:" << GenParticle[myT].pt << " - eta:" << GenParticle[myT].eta << endl
	       << "bbar.........: " << mybbar << "   pdg:" << GenParticle[mybbar].pdgId << " (" << pdg[abs(GenParticle[mybbar].pdgId)] << ")" << " - pt:" << GenParticle[mybbar].pt << " - eta:" << GenParticle[mybbar].eta << endl
	       << "tagging jet..: " << myj << "   pdg:" << GenParticle[myj].pdgId << " (" << pdg[abs(GenParticle[myj].pdgId)] << ")" << " - pt:" << GenParticle[myj].pt << " - eta:" << GenParticle[myj].eta << endl
	       << "Top decay-b..: " << mytb << "   pdg:" << GenParticle[mytb].pdgId << " (" << pdg[abs(GenParticle[mytb].pdgId)] << ")" << " - pt:" << GenParticle[mytb].pt << " - eta:" << GenParticle[mytb].eta << endl
	       << "Top decay-W..: " << mytW << "   pdg:" << GenParticle[mytW].pdgId << " (" << pdg[abs(GenParticle[mytW].pdgId)] << ")" << " - pt:" << GenParticle[mytW].pt << " - eta:" << GenParticle[mytW].eta << endl
	       << "W decay-lbar.: " << myWl << "   pdg:" << GenParticle[myWl].pdgId << " (" << pdg[abs(GenParticle[myWl].pdgId)] << ")" << " - pt:" << GenParticle[myWl].pt << " - eta:" << GenParticle[myWl].eta << endl
	       << "W decay-nu_l.: " << myWnu << "   pdg:" << GenParticle[myWnu].pdgId << " (" << pdg[abs(GenParticle[myWnu].pdgId)] << ")" << " - pt:" << GenParticle[myWnu].pt << " - eta:" << GenParticle[myWnu].eta << endl
	       << "Top last seen: " << lastT << "   pdg:" << GenParticle[lastT].pdgId << " (" << pdg[abs(GenParticle[lastT].pdgId)] << ")" << " - pt:" << GenParticle[lastT].pt << " - eta:" << GenParticle[lastT].eta << endl;
	  	  
	  // let's make sure all decay products have been identified (correctly)
	  assert(myT > 0);
	  assert(mybbar > 0);
	  assert(myj  > 0);
	  assert(mytb  > 0);
	  assert(mytW  > 0);
	  assert(myWl  > 0);
	  assert(myWnu > 0);	  	  
	  //     top        extra bbar    tagging d   top decay b  top decay W  W decay l+   W decay nu_l
	  
	  
	  
	  // ---------------------
	  // -- fiducial volume definition --
	  /*
	  ATLAS
	  ------------------------
	  Electrons 	pt> 25 GeV and |eta|<2.5
	  Muons     	pt> 25 GeV and |eta|<2.5
	  Jets      	pt> 30 GeV and |eta|<2.5
			pt> 35 GeV if  2.75 < |eta|< 3.5
	  Lepton, Jets  DR > 0.4
	  ETMISS   	Et> 30 GeV 
	  Transverse W mass  mTW> 50 GeV
	  Lepton, jet *variable*


	   */
	  // ---------------------
	  


	  
	  // ---------------------
	  // -- fill histograms --
	  // ---------------------	  

	}

  stream.close();
  ofile.close();
  //cout << "events: " << nevt << endl;
  return 0;
}
