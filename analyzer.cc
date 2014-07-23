//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Wed Jul  9 15:44:25 2014 by mkanalyzer.py
// Author:      Halil Gamsizkan
//-----------------------------------------------------------------------------
#include "analyzer.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCut.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "Riostream.h"
#include "TLatex.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line
int r;
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
TH1F h1("TopPT","TOP Pt;P\_{t} (GeV);Events/2GeV",150,0,300) ;
TH1F h2("TopETA","TOP eta; #eta;Events",60,0,6) ;
TH1F h3("BbarPT","\_${b}$ Pt;P\_{t} (GeV);Events/2GeV",100,0,200) ;
TH1F h4("BbarETA","\_${b}$ eta; #eta;Events",60,0,6) ;
TH1F h5("TagJetPT","Tagging Jet Pt;P\_{t} (GeV);Events/2GeV",100,0,200) ;
TH1F h6("TagJetETA","Tagging Jet eta; #eta;Events",60,0,6) ;
TH1F h7("bTopPT","B Pt;P\_{t} (GeV);Arbitrary Units",100,0,200) ;
TH1F h8("bTopETA","B Eta; #eta;Arbitrary Units",60,0,6) ;
TH1F h9("W-PT","W Pt;P\_{t} (GeV);Arbitrary Units",150,0,300) ;
TH1F h10("W-ETA","W-Boson eta; #eta;Arbitrary Units",60,0,6) ;
TH1F h11("lepPT","Lepton Pt;P\_{t} (GeV);Events/2GeV",100,0,200) ;
TH1F h12("lepETA","Lepton eta; #eta;Events",60,0,6) ;
TH1F h13("nuPT","Neutrino Pt;P\_{t} (GeV);Events/2GeV",100,0,200);
TH1F h14("nuETA","Neutrino eta; #eta;Events",60,0,6);
TH1F h15("lTopPT","Last Top Pt;P\_{t} (GeV);Events/2GeV",150,0,300);
TH1F h16("lTopETA","Last Top eta; #eta;Events",60,-6,6);
TH1F h17("reTm","Reconstructed Top Mass ; Mass (GeV);Events/0.5GeV",100,150,200);
TH1F h18("tTm","Reconstructed Top Mass from its Energy Momentum; Mass (GeV);Events/0.5GeV",100,150,200);
TH1F h19("Tm","Real Top Mass; Mass (GeV);Events/0.5GeV",100,150,200);
TH1F h20("WRm","Reconstructed W Mass ; Mass (GeV);Events/0.5GeV",160,40,120);
TH1F h21("rWm","Reconstructed W Mass ; Mass (GeV);Events/0.5GeV",160,40,120);
TH1F h22("Wm","Real W Mass; Mass (GeV);Events/0.5GeV",160,40,120);
TH1F h23("wldeltaeta","W-lepton Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h24("wldeltaphi","W-lepton Delta Phi; Delta Phi (rad);Events",50,0,5);
TH1F h26("bbartdeltaeta","BBar-Top Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h27("bbartdeltaphi","BBar-Top Phi; Delta Phi (rad);Events",50,0,5);
TH1F h28("tagtdeltaeta","Tag-Top Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h29("tagtdeltaphi","Tag-Top Delta Phi; Delta Phi (rad);Events",50,0,5);
TH1F h30("Wnuedeltaeta","W-Neutrino Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h31("Wnuedeltaphi","W-Neutrino Delta Phi; Delta Phi (rad);Events",50,0,5);
TH1F h32("ltWdeltaeta","Last Top - W Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h33("ltWdeltaphi","Last Top - W Delta Phi; Delta Phi (rad);Events",50,0,5);
TH1F h34("ltBdeltaeta","Last Top - B Delta Eta; Delta Eta (rad);Events",50,0,5);
TH1F h35("ltBdeltaphi","Last Top - B Delta Phi; Delta Phi (rad);Events",50,0,5);
TH1F h68("lepCutPt","lepCutPt; P\_{t};Events",250,0,500);
TH1F h69("lepCutEta","lepCutEta; #eta ;Events",160,-8,8);



TH1F h36("bbart","Bbar-Top Delta R; Delta R (rad);Events",50,0,5);
TH1F h37("tagt","Tagged JEt-Top Delta R; Delta R (rad);Events",50,0,5);
TH1F h38("wldeltar","W-Lepton Delta R; Delta R (rad);Events",50,0,5);
TH1F h39("wneudeltar","W-Neutrino Delta R; Delta R (rad);Events",50,0,5);
TH1F h40("ltwdeltar","Top-W Delta R; Delta R (rad);Events",50,0,5);
TH1F h41("ltbdeltar","Top-B Delta R; Delta R (rad);Events",50,0,5);

TH1F h42("tagjetdeltar"," Tagjet Delta R; Delta R ;Events",500,0,5);
TH1F h43("tagjetdeltar<0.1"," Tagjet Delta R < 0.1; Delta R ;Events",500,0,5);
TH1F h44("tagjetdeltar<0.3"," Tagjet Delta R < 0.3; Delta R ;Events",500,0,5);
TH1F h45("tagjetdeltar<1"," Tagjet Delta R < 1; Delta R ;Events",500,0,5);
TH1F h46("bjetdeltar"," B-jet Delta R; Delta R ;Events",500,0,5);
TH1F h47("bjetdeltar<0.1"," B-jet Delta R < 0.1; Delta R ;Events",500,0,5);
TH1F h48("bjetdeltar<0.3"," B-jetDelta R < 0.3; Delta R ;Events",500,0,5);
TH1F h49("bjetdeltar<1"," B-jet Delta R < 1; Delta R ;Events",500,0,5);

TH1F h50("tagjeta","Tag-jet Eta;  Eta ;Events",100,-5,5);
TH1F h51("tagjphi","Tag-jet  Phi; Phi (rad);Events",80,-4,4);
TH1F h52("bjeteta","B-jet Eta; Eta ;Events",100,-5,5);
TH1F h53("bjetphi","B-jet Phi;  Phi (rad);Events",80,-4,4);
TH1F h54("tagjetetar01","Tag-jet R<0.1 Eta; Eta ;Events",100,-5,5);
TH1F h55("tagjetphir01","Tag-jet R<0.1 Phi; Phi (rad);Events",80,-4,4);
TH1F h56("bjetetar01","B-jet R<0.1 Eta; Eta ;Events",100,-5,5);
TH1F h57("bjetphir01","B-jet R<0.1 Phi; Phi (rad);Events",80,-4,4);

TH1F h58("tagjetetar01","Tag-jet R>0.1 Eta; Eta ;Events",100,-8,8);
TH1F h59("tagjetetar03","Tag-jet R>0.3 Eta; Eta ;Events",100,-8,8);
TH1F h60("bjetetar01","B-jet R>0.1 Eta; Eta ;Events",100,-8,8);
TH1F h61("bjetetar03","B-jet R>0.3 Eta; Eta ;Events",100,-8,8);
TH1F h62("TagJetPT01","Tagging Jet Pt R>0.1 ;P\_{t} (GeV);Events/1GeV",200,0,200);
TH1F h63("TagJetPT03","Tagging Jet Pt R>0.3 ;P\_{t} (GeV);Events/1GeV",200,0,200);
TH1F h64("BJetPT01","B Jet Pt R>0.1;P\_{t} (GeV);Events/1GeV",200,0,200);
TH1F h65("BJetPT03","B Jet Pt R>0.3;P\_{t} (GeV);Events/1GeV",200,0,200);

TH1F h66("leppt","Lepton Pt ;P\_{t} (GeV);Events/1GeV",200,0,200);   

TH2F *Wlphiphi  = new TH2F("WlPhiphi","W Phi vs l Phi ; W-Phi ;l-Phi ",80,-4,4,80,-4,4);
TH2F *Wletaeta  = new TH2F("WlEtaeta","W Eta vs l Eta ; W-Eta ;l-Eta ",80,-8,8,80,-8,8);
TH2F *Wneuphiphi  = new TH2F("WnuPhiphi","W Phi vs Neutrino Eta ; W-Phi ;Neutrino-Phi ",80,-4,4,80,-4,4);
TH2F *Wneuetaeta  = new TH2F("WnuEtaeta","W Eta vs Neutrino Eta ; W-Eta ;Neutrino-Eta ",80,-8,8,80,-8,8);
TH2F *tWphiphi  = new TH2F("tWPhiphi","Top Phi vs W Phi ; Top-Phi ;W-Phi ",80,-4,4,80,-4,4);
TH2F *tWetaeta  = new TH2F("tWEtaeta","Top Eta vs W Eta ; Top-Eta ;W-Eta ",80,-8,8,80,-8,8);
TH2F *tBphiphi  = new TH2F("tbPhiphi","Top Phi vs b Phi ; Top-Phi ;b-Phi ",80,-4,4,80,-4,4);
TH2F *tBetaeta  = new TH2F("tbEtaeta","Top Eta vs b Eta ; Top-Eta ;b-Eta ",80,-8,8,80,-8,8);
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
//int tjetdeltaphi[];
int j=GenJet.size();
double tagjetdeltaphi[j];
double tagjetdeltaeta[j];
double tagjetdeltar[j];
double bjetdeltaphi[j];
double bjetdeltaeta[j];
double bjetdeltar[j];
int tagjet=0 ;
int bjet=0 ;
	for (int c=0; c<j ; c++) 
		{
		tagjetdeltaphi[c]=((GenParticle[myj].phi)-(GenJet[c].phi));
		tagjetdeltaphi[c]=abs(tagjetdeltaphi[c]);
			if (tagjetdeltaphi[c]>3.14)
				{
					tagjetdeltaphi[c]= 6.28-(tagjetdeltaphi[c]);
				}
		tagjetdeltaeta[c]=((GenParticle[myj].eta)-(GenJet[c].eta));				
		tagjetdeltaeta[c]=abs(tagjetdeltaeta[c])	;
 		tagjetdeltar[c]= TMath::Sqrt((tagjetdeltaphi[c]*tagjetdeltaphi[c])+(tagjetdeltaeta[c]*tagjetdeltaeta[c]));
		//cout<<"Jetr.....:"  << tagjetdeltar[c] <<" jetphi  : "<<tagjetdeltaphi[c]<<" myjphi : "<<GenParticle[myj].phi<<" jetphi : " << GenJet[c].phi << " jeteta    : "<<tagjetdeltaeta[c]<< endl;
		}
		double  mintag = tagjetdeltar[0];
			for (int k = 1 ; k < j ; k++ ) 
    			{
       				if ( tagjetdeltar[k] < mintag ) 
        				{	
           					mintag = tagjetdeltar[k];
						tagjet=k ;	
         				}
			}
        	h42. Fill(mintag);
	if (mintag<0.1 )
                {
                h43.Fill(mintag);
//		int o ;
  //              o++;
//		cout<<" Number of Events deltar<0.1  : "<<o<<endl;
        	h54.Fill(GenJet[tagjet].eta);
                h55.Fill(GenJet[tagjet].phi);
		}	
	if (mintag<0.3 )
                {
            	h44.Fill(mintag);
	//	int z ;
	//	z++;
	//	cout<<" Number of Events deltar<0.3  : "<<z<<endl;
		}	
	if (mintag<1 )
		{
		h45.Fill(mintag);
	//	cout<<"  mintag   :  "<<minimumtag<<"  jettagdeltar  : "<<tagjetdeltar[tagjet]<<endl;
	//	cout<<" myjphi : "<<GenParticle[myj].phi<<" jetphi : " << GenJet[tagjet].phi <<  "  phi : " <<tagjetdeltaphi[tagjet]<<endl;
	//	cout<<" myjeta : "<<GenParticle[myj].eta<<" jeteta : " << GenJet[tagjet].eta <<	 "  eta ; " <<tagjetdeltaeta[tagjet]<<endl ;
		}
	if(mintag >0.1)
	{
	h58.Fill(GenParticle[myj].eta);
	h62.Fill(GenParticle[myj].pt);
	}
	if(mintag >0.3)
        {
        h59.Fill(GenParticle[myj].eta);
        h63.Fill(GenParticle[myj].pt);
        }   
        for (int b=0; b<j ; b++)
                {
                bjetdeltaphi[b]=((GenParticle[mytb].phi)-(GenJet[b].phi));
                bjetdeltaphi[b]=abs(bjetdeltaphi[b]);
                        if (bjetdeltaphi[b]>3.14)
                                {
                                        bjetdeltaphi[b]= 6.28-(bjetdeltaphi[b]);
                                }
                bjetdeltaeta[b]=((GenParticle[mytb].eta)-(GenJet[b].eta));
                        if(bjetdeltaeta[b]<0)
                                {
                                        bjetdeltaeta[b]=-bjetdeltaeta[b];
                                }
                bjetdeltar[b]= TMath::Sqrt((bjetdeltaphi[b]*bjetdeltaphi[b])+(bjetdeltaeta[b]*bjetdeltaeta[b]));
                //cout<<"Jetr.....:"  << bjetdeltar[b] <<" bjetphi  : "<<bjetdeltaphi[b]<<" mytbphi : "<<GenParticle[mytb].phi<<" bjetphi : " << 
                }
                double  minb = bjetdeltar[0];
                        for (int h = 1 ; h < j ; h++ )
                        {
                                if ( bjetdeltar[h] < minb )
                                        {
                                                minb = bjetdeltar[h];
                                                bjet=h ;
                                        }
                        }
                h46. Fill(minb);
	if (minb<0.1 )
                {
                h47.Fill(minb);
          //      int g ;
            //    g++;
              //  cout<<" Number of Events deltar<0.1  : "<<g<<endl;
		h56.Fill(GenJet[bjet].eta);
		h57.Fill(GenJet[bjet].phi);
        	}
	if (minb<0.3 )
                {
                h48.Fill(minb);
          //      int f ;
            //    f++;
              //  cout<<" Number of Events deltar<0.3  : "<<f<<endl;
        	}
	if (minb<1 )
                {
                h49.Fill(minb);
        //      cout<<"  mintag   :  "<<minimumtag<<"  jettagdeltar  : "<<tagjetdeltar[tagjet]<<endl;
        //      cout<<" myjphi : "<<GenParticle[myj].phi<<" jetphi : " << GenJet[tagjet].phi <<  "  phi : " <<tagjetdeltaphi[tagjet]<<endl;
        //      cout<<" myjeta : "<<GenParticle[myj].eta<<" jeteta : " << GenJet[tagjet].eta <<  "  eta ; " <<tagjetdeltaeta[tagjet]<<endl ;
                }
	if(minb >0.1)
        {
        h60.Fill(GenParticle[mytb].eta);
        h64.Fill(GenParticle[mytb].pt);
        }
        if(minb >0.3)
        {
        h61.Fill(GenParticle[mytb].eta);
        h65.Fill(GenParticle[mytb].pt);
        }


		h50.Fill(GenJet[tagjet].eta);
                h51.Fill(GenJet[tagjet].phi);
		h52.Fill(GenJet[bjet].eta);
                h53.Fill(GenJet[bjet].phi);




/*	  cout << "Decay information:\n"
	       << "Top..........: " << myT << "   pdg:" << GenParticle[myT].pdgId << " (" << pdg[abs(GenParticle[myT].pdgId)] << ")" << " - pt:" << 
GenParticle[myT].pt << " - eta:" << GenParticle[myT].eta << endl
//	       << "bbar.........: " << mybbar << "   pdg:" << GenParticle[mybbar].pdgId << " (" << pdg[abs(GenParticle[mybbar].pdgId)] << ")" << " - 
pt:" << GenParticle[mybbar].pt << " - eta:" << GenParticle[mybbar].eta << endl
//	       << "tagging jet..: " << myj << "   pdg:" << GenParticle[myj].pdgId << " (" << pdg[abs(GenParticle[myj].pdgId)] << ")" << " - pt:" << 
GenParticle[myj].pt << " - eta:" << GenParticle[myj].eta << endl
//	       << "Top decay-b..: " << mytb << "   pdg:" << GenParticle[mytb].pdgId << " (" << pdg[abs(GenParticle[mytb].pdgId)] << ")" << " - pt:" << 
GenParticle[mytb].pt << " - eta:" << GenParticle[mytb].eta << endl
	       << "Top decay-W..: " << mytW << "   pdg:" << GenParticle[mytW].pdgId << " (" << pdg[abs(GenParticle[mytW].pdgId)] << ")" << " - pt:" << 
GenParticle[mytW].pt << " - eta:" << GenParticle[mytW].eta << endl
	       << "W decay-lbar.: " << myWl << "   pdg:" << GenParticle[myWl].pdgId << " (" << pdg[abs(GenParticle[myWl].pdgId)] << ")" << " - pt:" << GenParticle[myWl].pt << " - eta:" << GenParticle[myWl].eta << endl
	       << "W decay-nu_l.: " << myWnu << "   pdg:" << GenParticle[myWnu].pdgId << " (" << pdg[abs(GenParticle[myWnu].pdgId)] << ")" << " - pt:" << GenParticle[myWnu].pt << " - eta:" << GenParticle[myWnu].eta << endl
		<< "Top last seen: " << lastT << "   pdg:" << GenParticle[lastT].pdgId << " (" << pdg[abs(GenParticle[lastT].pdgId)] << ")" << " - pt:" << GenParticle[lastT].pt << " - eta:" << GenParticle[lastT].eta << endl;
*/	  	  
	  // let's make sure all decay products have been identified (correctly)
	  assert(myT > 0);
	  assert(mybbar > 0);
	  assert(myj  > 0);
	  assert(mytb  > 0);
	  assert(mytW  > 0);
	  assert(myWl  > 0);
	  assert(myWnu > 0);	  	  
	  //     top        extra bbar    tagging d   top decay b  top decay W  W decay l+   W decay nu_l
	
 		//Top Mass Reconstruction
	float bpx = GenParticle[mytb].px ;
	float bpy = GenParticle[mytb].py ;
	float bpz = GenParticle[mytb].pz ;
        float bE = GenParticle[mytb].energy ;
        float Wpx = GenParticle[mytW].px ;
        float Wpy = GenParticle[mytW].py ;
        float Wpz = GenParticle[mytW].pz ;   
        float WE = GenParticle[mytW].energy ;
        float TE = GenParticle[lastT].energy ;   
        float Tpx = GenParticle[lastT].px ;
        float Tpy = GenParticle[lastT].py ;                               
        float Tpz = GenParticle[lastT].pz ;                               
        float bwpx = (bpx+Wpx);
	float bwpy = (bpy+Wpy);
        float bwpz = (bpz+Wpz);
        float bwE = (bE+WE);
		//W Mass Reconstruction
	float nupx = GenParticle[myWnu].px ;
        float nupy = GenParticle[myWnu].py ;
        float nupz = GenParticle[myWnu].pz ;
        float nuE = GenParticle[myWnu].energy ;
        float lE = GenParticle[myWl].energy ;
        float lpx = GenParticle[myWl].px ;
        float lpy = GenParticle[myWl].py ;
        float lpz = GenParticle[myWl].pz ;
        float nulpx = (nupx+lpx);
        float nulpy = (nupy+lpy);
        float nulpz = (nupz+lpz);
        float nulE = (nuE+lE);


//W-l deltaeta
        float Wldeltaeta = ((GenParticle[mytW].eta)-(GenParticle[myWl].eta));
        if(Wldeltaeta < 0){Wldeltaeta = -Wldeltaeta;}
//W-l deltaphi
        float Wldeltaphi ;
        
                if((GenParticle[mytW].phi)>(GenParticle[myWl].phi))
                {
                         Wldeltaphi = ((GenParticle[mytW].phi)-(GenParticle[myWl].phi));
                }
                if((GenParticle[myWl].phi)>(GenParticle[mytW].phi))
                {
                Wldeltaphi=((GenParticle[myWl].phi)-(GenParticle[mytW].phi));
                }
                if((Wldeltaphi>3.14))
                {
                Wldeltaphi=(6.28)-Wldeltaphi ;
                }
//bbar-top deltaeta
        float bbartdeltaeta = ((GenParticle[mybbar].eta)-(GenParticle[myT].eta));
        if(bbartdeltaeta < 0){bbartdeltaeta = -bbartdeltaeta;}
//bbar-top deltaphi
        float bbartdeltaphi ;
        
                if((GenParticle[myT].phi)>(GenParticle[mybbar].phi))
                {
                         bbartdeltaphi = ((GenParticle[myT].phi)-(GenParticle[mybbar].phi));
                }
                if((GenParticle[mybbar].phi)>(GenParticle[myT].phi))
                {
                bbartdeltaphi=((GenParticle[mybbar].phi)-(GenParticle[myT].phi));
                }
                if((bbartdeltaphi>3.14))
                {
                bbartdeltaphi=(6.28)-bbartdeltaphi ;
                }
//tagj-top deltaeta
        float tagtdeltaeta = ((GenParticle[myj].eta)-(GenParticle[myT].eta));
        if(tagtdeltaeta < 0){tagtdeltaeta = -tagtdeltaeta;}
//tagj-top deltaphi
        float tagtdeltaphi ;
        
                if((GenParticle[myT].phi)>(GenParticle[myj].phi))
                {
                         tagtdeltaphi = ((GenParticle[myT].phi)-(GenParticle[myj].phi));
                }
                if((GenParticle[myj].phi)>(GenParticle[myT].phi))
                {
                tagtdeltaphi=((GenParticle[myj].phi)-(GenParticle[myT].phi));
                }
                if((tagtdeltaphi>3.14))
                {
                tagtdeltaphi=(6.28)-tagtdeltaphi ;
                }
//W-neu deltaeta
        float Wneudeltaeta = ((GenParticle[mytW].eta)-(GenParticle[myWnu].eta));
        if(Wneudeltaeta < 0){Wneudeltaeta = -Wneudeltaeta;}
//W-neu deltaphi
        float Wneudeltaphi ;
        
                if((GenParticle[mytW].phi)>(GenParticle[myWnu].phi))
                {
                         Wneudeltaphi = ((GenParticle[mytW].phi)-(GenParticle[myWnu].phi));
                }
                if((GenParticle[myWnu].phi)>(GenParticle[mytW].phi))
                {
                Wneudeltaphi=((GenParticle[myWnu].phi)-(GenParticle[mytW].phi));
                }
                if((Wneudeltaphi>3.14))
                {
                Wneudeltaphi=(6.28)-Wneudeltaphi ;
                }
//lastTop-W deltaeta
        float ltWdeltaeta = ((GenParticle[mytW].eta)-(GenParticle[lastT].eta));
        if(ltWdeltaeta < 0){ltWdeltaeta = -ltWdeltaeta;}
//lastTop-W deltaphi 
        float ltWdeltaphi ;
        
                if((GenParticle[mytW].phi)>(GenParticle[lastT].phi))
                {  
                         ltWdeltaphi = ((GenParticle[mytW].phi)-(GenParticle[lastT].phi));
                }
                if((GenParticle[lastT].phi)>(GenParticle[mytW].phi))
                {
                ltWdeltaphi=((GenParticle[lastT].phi)-(GenParticle[mytW].phi));
                }
                if((ltWdeltaphi>3.14))
                {
                ltWdeltaphi=(6.28)-ltWdeltaphi ;
                }
//lastTop-B deltaeta
        float ltBdeltaeta = ((GenParticle[mytb].eta)-(GenParticle[lastT].eta));
        if(ltBdeltaeta < 0){ltBdeltaeta = -ltBdeltaeta;}
//lastTop-B deltaphi
        float ltBdeltaphi ;
 
                if((GenParticle[mytb].phi)>(GenParticle[lastT].phi))
                {
                         ltBdeltaphi = ((GenParticle[mytb].phi)-(GenParticle[lastT].phi));
                }
                if((GenParticle[lastT].phi)>(GenParticle[mytb].phi))
                {
                ltBdeltaphi=((GenParticle[lastT].phi)-(GenParticle[mytb].phi));
                }
                if((ltBdeltaphi>3.14))
                {
                ltBdeltaphi=(6.28)-ltBdeltaphi ;
                }
 //     T MASS and W MASS Reconstruction
	float recTmass = TMath::Sqrt((bwE*bwE)-(bwpx*bwpx)-(bwpy*bwpy)-(bwpz*bwpz));
	float tTmass = TMath::Sqrt((TE*TE)-(Tpx*Tpx)-(Tpy*Tpy)-(Tpz*Tpz)) ;
//	  cout << "MASSTop..........: " << GenParticle[myT].mass  << "recTmass....: " << recTmass << "tTmass......: " << tTmass <<endl ;
//	  cout << "Deltaeta ... : " << Wldeltaphi <<endl ;
	float Wldeltar= TMath::Sqrt((Wldeltaphi*Wldeltaphi)+(Wldeltaeta*Wldeltaeta));
	float bbartdeltar=TMath::Sqrt((bbartdeltaphi*bbartdeltaphi)+(bbartdeltaeta*bbartdeltaeta));
	float tagtdeltar= TMath::Sqrt((tagtdeltaphi*tagtdeltaphi)+(tagtdeltaeta*tagtdeltaeta));
	float Wneudeltar= TMath::Sqrt((Wneudeltaphi*Wneudeltaphi)+(Wneudeltaeta*Wneudeltaeta));
	float ltWdeltar= TMath::Sqrt((ltWdeltaphi*ltWdeltaphi)+(ltWdeltaeta*ltWdeltaeta));
	float ltBdeltar= TMath::Sqrt((ltBdeltaphi*ltBdeltaphi)+(ltBdeltaeta*Wldeltaeta));
	float recWmass = TMath::Sqrt((nulE*nulE)-(nulpx*nulpx)-(nulpy*nulpy)-(nulpz*nulpz));
        float tWmass = TMath::Sqrt((WE*WE)-(Wpx*Wpx)-(Wpy*Wpy)-(Wpz*Wpz)) ;
//	cout << "MASSW..........: " << GenParticle[mytW].mass  << "recWmass....: " << recWmass << "tWmass......: " << tWmass <<endl ;
	
//	 eta[mytb]=GenParticle[mytb].eta ;
//	 phi[mytb]=GenParticle[mytb].phi ;
//2D PHI-PHI
	Wlphiphi->Fill(GenParticle[mytW].phi,GenParticle[myWl].phi);
	Wletaeta->Fill(GenParticle[mytW].eta,GenParticle[myWl].eta);
	Wneuphiphi->Fill(GenParticle[mytW].phi,GenParticle[myWnu].phi);
        Wneuetaeta->Fill(GenParticle[mytW].eta,GenParticle[myWnu].eta);
	tWphiphi->Fill(GenParticle[lastT].phi,GenParticle[mytW].phi);
        tWetaeta->Fill(GenParticle[lastT].eta,GenParticle[mytW].eta);
	tBphiphi->Fill(GenParticle[lastT].phi,GenParticle[mytb].phi);
        tBetaeta->Fill(GenParticle[lastT].eta,GenParticle[mytb].eta);
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
//int r ;
	if( GenParticle[myWl].pt>20 && (abs(GenParticle[myWl].eta))<2.5 
&& GenParticle[myj].pt>20 && GenParticle[mytb].pt>20 && 
GenMET[0].pt>20 && (abs(GenParticle[myj].eta))<4 && 
(abs(GenParticle[mytb].eta))<4)
	{
	h66. Fill(GenParticle[myWl].pt);
//	int r ;
	r++ ;
//	double k=(r/nevents);
//	cout<<k<<"  sadad   "<<r<<endl;
	} 
	if(GenParticle[myWl].pt>20 && (abs(GenParticle[myWl].eta)<2.5))
	{
	h68. Fill(GenParticle[myWl].pt);
	h69. Fill(abs(GenParticle[myWl].eta));
	} 

	h1. Fill(GenParticle[myT].pt);
	h2. Fill(abs(GenParticle[myT].eta));
	 
        h3. Fill(GenParticle[mybbar].pt);
        h4. Fill(abs(GenParticle[mybbar].eta));
	 
        h5. Fill(GenParticle[myj].pt);
        h6. Fill(abs(GenParticle[myj].eta));
	 
        h7. Fill(GenParticle[mytb].pt);
        h8. Fill(abs(GenParticle[mytb].eta));
	 
        h9. Fill(GenParticle[mytW].pt);
        h10. Fill(abs(GenParticle[mytW].eta));
	 
        h11. Fill(GenParticle[myWl].pt);
        h12. Fill(abs(GenParticle[myWl].eta));
	 
        h13. Fill(GenParticle[myWnu].pt);
        h14. Fill(abs(GenParticle[myWnu].eta));
	 
        h15. Fill(GenParticle[lastT].pt);
        h16. Fill(abs(GenParticle[lastT].eta));
	
	h17. Fill(recTmass);
        h18. Fill(tTmass);
        h19. Fill(GenParticle[myT].mass);  
	h20. Fill(recWmass);
        h21. Fill(tWmass);
        h22. Fill(GenParticle[mytW].mass);
	h23. Fill(Wldeltaeta);
	h24. Fill(Wldeltaphi);
//	h25. Fill(Wldeltar);
	h26. Fill(bbartdeltaeta);
        h27. Fill(bbartdeltaphi);
	h28. Fill(tagtdeltaeta);
        h29. Fill(tagtdeltaphi);
	h30. Fill(Wneudeltaeta);
        h31. Fill(Wneudeltaphi);
	h32. Fill(ltWdeltaeta);
        h33. Fill(ltWdeltaphi);
	h34. Fill(ltBdeltaeta);
        h35. Fill(ltBdeltaphi);
h36. Fill(bbartdeltar);
h37. Fill(tagtdeltar);
h38. Fill(Wldeltar);
h39. Fill(Wneudeltar);
h40. Fill(ltWdeltar);
h41. Fill(ltBdeltar);	
//h42. Fill(minimum);
//hpxpy->SetLineColor(kRed);
//h43. Fill(minimumtag);	
  // ---------------------
	  // -- fill histograms --
	  // ---------------------	  

}
float k=r/17000.00;
cout<<k<<" sdadsad   "<<r<<endl;
  stream.close();
  ofile.close();
  //cout << "events: " << nevt << endl;
  return 0;
}

