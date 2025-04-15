//TOF
#include <climits>

#include "Pythia8/Pythia.h"
using namespace Pythia8; 

#include "qbh.h"
using namespace QBH; 

/*****************************************************************************/
/* Event generator code (main program).                                      */
/*****************************************************************************/
int main(int argc, char *argv[])
{
  bool status;

  // Output LHEF file.
  const bool OUTPUT = true;

  // Declare Pythia object.
  Pythia* pythia = new Pythia();

  // Declare and initialize quantum black hole object.
  bool initialize = true;
  QuantumBlackHole* qbh = new QuantumBlackHole(pythia,initialize);

  // Decalare derived LHA user process object.
  dLHAup* lhaPtr = new dLHAup;

  pythia->setLHAupPtr(lhaPtr);  
  pythia->readString("Beams:frameType = 5");

  string LHEFilename;
  string QBHFilename;
  LHEFilename = "LHEFfile.lhe";
  QBHFilename = "QBHfile.lhe";

  // Number of event, black hole charge, and inital state.
  if (argc == 2) {
    pythia->readString("Main:NumberOfEvents = 20000");
    qbh->setMplanck(atof(argv[1]));
    qbh->setQstate(4);
    qbh->setIstate(0);
  }
  else if (argc == 3) {
    pythia->readString("Main:NumberOfEvents = 0");
    qbh->setQstate(atoi(argv[1]));
    qbh->setIstate(atoi(argv[2]));
  }
  else if (argc == 4) {
    pythia->readString("Main:NumberOfEvents = " + (string)argv[1]);
    qbh->setQstate(atoi(argv[2]));
    qbh->setIstate(atoi(argv[3]));
  }
  else if (argc ==6) {
    pythia->readString("Main:NumberOfEvents = " + (string)argv[1]);
    qbh->setQscale(false);
    qbh->setYRform(false);
    qbh->setTrap(false);
    qbh->setSM(true);
    qbh->setMajorana(false);
    qbh->setChiral(false);
    qbh->setHiggs(false);
    qbh->setGraviton(false);
    qbh->setTotdim(4+atoi(argv[2]));
    qbh->setPlanckdef(1);
    qbh->setRS1(true);
    qbh->setMplanck(atoi(argv[3]));
    qbh->setIstate(3);
    qbh->setQstate(9);
    qbh->setIstate(3);
    qbh->setMinmass(atoi(argv[4]));
    qbh->setMaxmass(13000);
    qbh->setEcm(13000);
    LHEFilename = "/afs/cern.ch/work/z/zhangj/private/CMSSW_7_4_15/src/QBH/LHEFfile"+(string)argv[5]+".lhe";
    QBHFilename = "/afs/cern.ch/work/z/zhangj/private/CMSSW_7_4_15/src/QBH/QBHfile"+(string)argv[5]+".lhe";
  }
  else {
    qbh->setQstate(9);
    qbh->setIstate(3);
  }

  
  // Initialize PDFs.
  qbh->setLHAglue(10042);  //CTEQ6L1
  // qbh->setLHAglue(192800); //NNPDF21_100
  // qbh->setLHAglue(10800);  //CT10
  // qbh->setLHAglue(10550);  //CTEQ6.6
  // qbh->setLHAglue(19050);  //CTEQ5M
  // qbh->setLHAglue(21000);  //MSTW2008lo
  // qbh->setLHAglue(29041);  //MRST98lo fit
  // qbh->setLHAglue(0);      //PYTHIA internal

  // Set some PYTHIA switches.
  (void)pythia->readString("HardQCD:all = off");
  (void)pythia->readString("SoftQCD:nonDiffractive = off");
  (void)pythia->readString("SoftQCD:all       = off");
  (void)pythia->readString("PartonLevel:ISR   = on");
  (void)pythia->readString("PartonLevel:FSR   = on");
  (void)pythia->readString("PartonLevel:MPI    = on");
  (void)pythia->readString("HadronLevel:all   = off");
  (void)pythia->readString("Check:history     = on");
  (void)pythia->readString("Check:nErrList    = 10");

  if (OUTPUT) (void)lhaPtr->openLHEF(LHEFilename);

  // Initialize Pythia object.
  status = pythia->init();

  if (OUTPUT) (void)lhaPtr->initLHEF();  
  delete qbh;

  // Loop over events.
  if (status) {
    for (int iEvent=0;iEvent<pythia->mode("Main:numberOfEvents"); ++iEvent) {

      // Generate event.
      status = pythia->next();
      if (!status) continue;

      // Event listing.
      if (iEvent < pythia->mode("Next:numberShowEvent")) {
        //pythia->info.list();
        pythia->process.list();
        //pythia->event.list();
      }

      // Process analysis.
      //for (int i=0;i<pythia->process.size();i++) {
        // Black hole.
        //if (pythia->process[i].id() == 40) {
        //  double px = pythia->process[3].px() + pythia->process[4].px();
        //  double py = pythia->process[3].py() + pythia->process[4].py();
        //  double pz = pythia->process[3].pz() + pythia->process[4].pz();
        //  double e  = pythia->process[3].e()  + pythia->process[4].e();
        //  fprintf(stream,"%f %f %f %f %f\n",px,py,pz,e,pythia->process[i].m());
        //}
      //}
      // PDG ID codes.
      //fprintf(stream,"%d %d\n",pythia->process[5].id(),pythia->process[6].id());
      // Particle kinematics.
      //fprintf(stream,"%f %f %f %f %f\n",
      //               pythia->process[5].px(),
      //               pythia->process[5].py(),
      //               pythia->process[5].pz(),
      //               pythia->process[5].e(),
      //               pythia->process[5].m());
      //   fprintf(stream,"%f %f %f %f %f\n",
      //               pythia->process[6].px(),
      //               pythia->process[6].py(),
      //               pythia->process[6].pz(),
      //               pythia->process[6].e(),
      //               pythia->process[6].m());

      if (OUTPUT) (void)lhaPtr->eventLHEF(false);
    } // End of event loop.
  }
   
  // Termination.
  pythia->stat();
  
  if (OUTPUT) {
    (void)lhaPtr->closeLHEF();
    // Write QBH banner information to QBHfile.lhe
    qbh->trailer(LHEFilename, QBHFilename);
  }

  // Time job exectution.
  long t = clock();
  double dt = (double)t;
  if (t < 0) dt += 2.0 * (double)LONG_MAX;
  long min = (long)((dt/(double)CLOCKS_PER_SEC) / 60.0);
  long sec = (long)(dt/(double)CLOCKS_PER_SEC) % 60;
  printf("\nProcessing time %ld min %ld sec\n",min,sec);

  if (status) {
    printf("Normal sucessful completion.\n\n");
  }
  else {
    printf("Possible execution problem.\n\n");
  }
  return status;
}
/*****************************************************************************/
//EOF

