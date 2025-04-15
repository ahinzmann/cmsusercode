//TOF
/*****************************************************************************/
/* QBH version 3.02 - Quantum Black Hole Generator                           */
/* D.M. Gingrich 30-April-2022                                               */
/*****************************************************************************/
#include "Pythia8/Pythia.h"
using namespace Pythia8; 

#include "qbh.h"

namespace QBH {

/*****************************************************************************/
/* Decalare static variables to make them global.                            */
/* Give reasonable inital values.                                            */
/*****************************************************************************/
bool   QuantumBlackHole::saveQscale    = true;
bool   QuantumBlackHole::savePoisson   = false;
bool   QuantumBlackHole::saveYRform    = false;
bool   QuantumBlackHole::saveTrap      = false;
bool   QuantumBlackHole::saveSM        = false;
bool   QuantumBlackHole::saveHiggs     = true;
bool   QuantumBlackHole::saveGraviton  = true;
bool   QuantumBlackHole::saveMajorana  = false;
bool   QuantumBlackHole::saveChiral    = false;
bool   QuantumBlackHole::saveRS1       = false;
int    QuantumBlackHole::saveTotdim    = 10;
int    QuantumBlackHole::savePlanckdef = 3;
int    QuantumBlackHole::saveQstate    = 9;
int    QuantumBlackHole::saveIstate    = 3;
int    QuantumBlackHole::saveLHAglue   = 10042;
double QuantumBlackHole::saveEcm       = 13000.0;
double QuantumBlackHole::saveMplanck   =  3000.0;
double QuantumBlackHole::saveMinmass   = saveMplanck;
double QuantumBlackHole::saveMaxmass   = 3.0*saveMplanck;
int    QuantumBlackHole::fsize[]={91,44,41,42,41,47,41,48};
double QuantumBlackHole::br40[2];
double QuantumBlackHole::br30[7];
double QuantumBlackHole::br21[6];
double QuantumBlackHole::br20[3];
double QuantumBlackHole::br11[6];
double QuantumBlackHole::br10[3];
double QuantumBlackHole::br00[18];
double QuantumBlackHole::br02[14];
double QuantumBlackHole::za[91];
double QuantumBlackHole::ya[91];
PDF*   QuantumBlackHole::pdfPtr = NULL;
Pythia* QuantumBlackHole::pdg = NULL;
Rndm*   Random::random = NULL;

/*****************************************************************************/
/* Class Random member function delerations.                                 */
/*****************************************************************************/
Random::Random(void)
{
  random = new Rndm();
  random->init(363685);
}
/*****************************************************************************/
/* Random number generator greater than 0.0 to less than 1.0                 */
/*****************************************************************************/
inline double Random::ran(void)
{
  return random->flat();
}
/*****************************************************************************/
/* Returns true with probability a.                                          */
/*****************************************************************************/
inline bool Random::logic(const double a) 
{
  if (Random::ran() > a) {
    return false;
  }
  else {
    return true;
  }
}
/*****************************************************************************/
/* Uniform random number in range [a,b]                                      */
/*****************************************************************************/
inline double Random::uniform(const double a,const double b)
{
  return (a + Random::ran()*(b-a));
}
/*****************************************************************************/
/* Random up-quark PDG code: u = 2, c = 4, t = 6.                            */
/*****************************************************************************/
inline int Random::ranup(void)
{
  return (2*(int)(3*Random::ran()) + 2);
}
/*****************************************************************************/
/* Random down-quark PDG code: d = 1, s = 3, b = 5.                          */
/*****************************************************************************/
inline int Random::randn(void)
{
  return (2*(int)(3*Random::ran()) + 1);
}
/*****************************************************************************/
/* Random negative charged lepton PDG code: e = 11, mu = 13, tau = 15.       */
/*****************************************************************************/
inline int Random::ranlp(void)
{
  return (2*(int)(3*Random::ran()) + 11);
}
/*****************************************************************************/
/* Random neutrino PDG code: nu_e = 12, nu_mu = 14, nu_tua = 16,             */
/*                      nu-bar_e = -12, nu-bar_mu = -14, nu-bar_tau = -16.   */
/*****************************************************************************/
inline int Random::rannu(void)
{
  int id = 2*(int)(3*Random::ran()) + 12;
  return id;
}
/*****************************************************************************/
/* Class Decay member function delerations.                                  */
/*****************************************************************************/
/* Randomly rotated 2-vector (px,py) of length pt.                           */
/*****************************************************************************/
void Decay::randomAzimuth(const double pt,double& px,double& py)
{
  double c, s, cs;

  do {
    c = 2.0*Random::ran() - 1.0;
    s = 2.0*Random::ran() - 1.0;
    cs = c*c + s*s;
  }
  while ((cs > 1.0) || (cs == 0.0));

  double qt = pt / cs;
  px = (c*c - s*s) * qt;
  py = 2.0 * c * s * qt;

  return;
}
/*****************************************************************************/
/* Transforms pp (given in rest frame of ps) into pf (in lab).               */
/* N.B. p(0,1,2,3,4) = (px,py,pz,e,m)                                        */
/*****************************************************************************/
void Decay::boost(const double* ps,double* pf)
{
  if (ps[3] == ps[4]) {
    for (int i=0;i<4;i++) {
      pf[i] = pp[i];
    }
  }
  else {
    double pf4 = (pp[0]*ps[0] + pp[1]*ps[1] + pp[2]*ps[2] + pp[3]*ps[3]) 
               / ps[4];
    double fn  = (pf4 + pp[3]) / (ps[3] + ps[4]);

    for (int i=0;i<3;i++) {
      pf[i] = pp[i] + fn*ps[i];
    }
    if (pp[4] > 1.0e-5) {
      pf[3] = pf4;
    }
    else {
      pf[3] = sqrt(pf[0]*pf[0]+pf[1]*pf[1]+pf[2]*pf[2]);
    }
  }

  pf[4] = pp[4];
}
/*****************************************************************************/
/* r is rotation matrix to get from vector p to z-axis, followed by a        */
/* rotation by psi about z-axis, where cp = cos(-psi), sp = sin(-psi).       */
/*****************************************************************************/
void Decay::matrix(const double* p,const double cp,const double sp)
{
  double ct, st, cf, sf;
  const double ptcut = 1.0e-20;
  const double wn = 1.0;

  double pt = p[0]*p[0] + p[1]*p[1];
  double pm = p[2]*p[2] + pt;

  if (pt <= pm*ptcut) {
    ct = wn * p[2]/fabs(p[2]);
    st = 0.0;
    cf = 1.0;
    sf = 0.0;
  }
  else {
    pm = sqrt(pm);
    pt = sqrt(pt);
    ct = p[2] / pm;
    st = pt   / pm;
    cf = p[0] / pt;
    sf = p[1] / pt;
  }

  r[0][0] =  cp*cf*ct + sp*sf;
  r[0][1] =  cp*sf*ct - sp*cf;
  r[0][2] = -cp*st;
  r[1][0] = -cp*sf + sp*cf*ct;
  r[1][1] =  cp*cf + sp*sf*ct;
  r[1][2] = -sp*st;
  r[2][0] =  cf*st;
  r[2][1] =  sf*st;
  r[2][2] =  ct;

  return;
}
/*****************************************************************************/
/* Rotates a 3-vector by the inverse of rotation matrix r.                   */
/*****************************************************************************/
void Decay::rotate(void)
{
  double s0 = pp[0]*r[0][0] + pp[1]*r[1][0] + pp[2]*r[2][0];
  double s1 = pp[0]*r[0][1] + pp[1]*r[1][1] + pp[2]*r[2][1];
  double s2 = pp[0]*r[0][2] + pp[1]*r[1][2] + pp[2]*r[2][2];
  pp[0] = s0;
  pp[1] = s1;
  pp[2] = s2;

  return;
}
/*****************************************************************************/
/* Generates decay 0 -> 1 + 2                                                */
/* pcm is C.M. momentum or final-state particles.                            */
/* costh = cos(theta) in p0 rest frame (>1 for isotropic).                   */
/*****************************************************************************/
void Decay::twoBody(const double* p0,double* p1,double* p2,
                    const double pcm,const bool zaxis)
{
  // Choose C.M. angles.
  double cos = Random::uniform(-1.0,1.0);
  double sin = sqrt(1.0 - cos*cos);
  randomAzimuth(pcm*sin,pp[0],pp[1]);

  // pp is momentum of 2 in C.M.
  pp[2] = -pcm * cos;
  pp[3] = sqrt(p2[4]*p2[4] + pcm*pcm);
  pp[4] = p2[4];

  // Rotate if necessary.
  if (!zaxis) {
    matrix(p0,1.0,0.0);
    rotate();
  }

  // Boost from C.M. to lab frame.
  boost(p0,p2);
  for (int i=0;i<4;i++) {
    p1[i] = p0[i] - p2[i];
  }

  return;
}
/*****************************************************************************/
/* Class BlackHole member function delerations.                              */
/*****************************************************************************/
/* Initialize the parton distibutions.                                       */
/*****************************************************************************/
void QuantumBlackHole::setLHAglue(int gluecode)
{
  Info info;

  // Decalare PDF object to initalize PDFs.
  saveLHAglue = gluecode;

  // LHAPDF external.
  if      (gluecode == 10042) {
    string pdfSet = "LHAPDF6:cteq6l1";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if      (gluecode == 10550) {
    string pdfSet = "LHAPDF6:cteq66.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 19050) {
    string pdfSet = "LHAPDF6:cteq5m.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 10800) {
    string pdfSet = "LHAPDF6:CT10.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 29041) {
    string pdfSet = "LHAPDF6:MRST98lo.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 21000) {
    string pdfSet = "LHAPDF6:MSTW2008lo68cl.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 192800) {
    string pdfSet = "LHAPDF6:NNPDF21_100.LHgrid";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 303600) {
    string pdfSet = "LHAPDF6:NNPDF31_nnlo_as_0118";
    pdfPtr = new LHAPDF(2212,pdfSet,&info);
  }
  else if (gluecode == 0) {
    // PYTHIA internal.
    pdfPtr = new CTEQ5L(2212);
  }
  else {
    (void)printf("No PDF selected\n.");
    exit(EXIT_FAILURE);
  }
}
/*****************************************************************************/
/* Read in file for calculation of trapped surface cross section.            */
/*****************************************************************************/
void QuantumBlackHole::setTrap(bool inelast)
{
  const bool output = false;
  saveTrap = inelast;

  // If radiation allowed read in lower bounds.
  if (inelast) {

    string filename;
    ostringstream convert;
    convert << totdim();
    filename = "./MLB_D" + convert.str() + ".data";
    cout << filename;
    ifstream inFile;
    inFile.open(filename.c_str());
    if (!inFile) {
      cerr << "QBH::setTrap: Unable to open file" << endl;
    }

    const int ASIZE = fsize[totdim()-4];
    for (int i=0;i<ASIZE;i++) {
      inFile >> za[i] >> ya[i]; 
    }
    for (int i=0;i<ASIZE;i++) {
      za[i] /= za[ASIZE-1];
      if (output) {
        printf("%f %f\n",za[i],ya[i]); 
      }
    }

    inFile.close();
  }
}
/*****************************************************************************/
/* Give all data members of the black hole class sensible values.            */
/*****************************************************************************/
QuantumBlackHole::QuantumBlackHole(Pythia* pythia,bool initialize)
{
  if (!initialize) return;

  // Save Pythia point to access particle data table later.
  pdg = pythia;

  // Construct ramdom number generator.
  Random();

  // Proton-proton CMS energy.
  setEcm(13000.0);

  // Planck scale.
  setMplanck(3000.0);

  // Minimum black hole mass.
  setMinmass(mplanck());

  // Maximum black hole mass.
  const double ZETA=3.0;
  if (ZETA*mplanck() < ecm()) {
    setMaxmass(ZETA*mplanck());
  }
  else {
    setMaxmass(ecm());
  }

  // Total number of dimensions.
  if (RS1()) {
    setTotdim(5);
  }
  else {
    setTotdim(10);
  }

  // Definition of Planck scale.
  if (RS1()) {
    // Randall and Sundrum definition of Planck scale.
    setPlanckdef(1);
  }
  else {
    // Dimopoulos and Landsberg definition of Planck scale.
    // setPlanckdef(2);
    // PDG definition of Planck scale.
    setPlanckdef(3);
  }

  // Definition of QCD scale for parton distibutions.
  // QCD scale is mass in PDFs.
  // setQscale(false);
  // QCD scale is 1/radius in PDFs.
  setQscale(true);
  
  // Include two-particle Poisson probabilities in cross section.
  // setPoisson(true);
  setPoisson(false);

  // Yoshino-Rychkov form factor.
  // setYRform(true);
  setYRform(false);

  // Which particles to include in decays.
  // Included Higgs.
  setHiggs(true);
  // Include gravition.
  setGraviton(true);
  // Discovered SM particles only.
  //setGraviton(false);

  // Impose Standard Model global symmetries.
  //setSM(true);
  // Allow violation of global symmetries.
  setSM(false);

  // Neutrinos are both left- and right-handed.
  // setChiral(true);
  // Neutrinos are only left-handed.
  setChiral(false);

  // Majorana or Dirac neutrinos.
  // Majorana neutrinos.
  // setMajorana(true);
  // Dirac neutrinos.
  setMajorana(false);

  // Quantum black hole state variables.
  setQstate(9);
  setIstate(3);

  // Totally inelastic cross section.
  setTrap(false);
  // Trapped surface lower-bound cross section.
  // setTrap(false);

  // Calculate branching ratios.
  const bool output = false;
  double ql, nl, qn, ll, nn, WH, WG, qH, qG, gH, gG, HH, pH, ZH, pG, ZG, GG, HG;

  if (sm()) {
    ql = 0;
    qn = 0;
    nl = 3;
    ll = 3;
    if (chiral()) {
      nn = 3;
    }
    else {
      nn = 3 * 3./4.;
    }
  }
  else {
    ql =  9;
    ll =  9;
    if (majorana()) {
      qn = 9;
      nl = 9;
      if (chiral()) {
        nn = 6;
      }
      else {
        nn = 6 * 1./4.;
      }
    }
    else {
      qn = 18;
      nl = 18;
      if (chiral()) {
        nn = 21;
      }
      else {
        nn = 21 * 3./4.;
      }
    }
  }
  if (higgs()) {
    WH = 1 * 3./4.;
    qH = 3 * 1./3.;
    gH = 1 * 3./4.;
    HH = 1 * 1./4.;
    pH = 1 * 3./4.;
    ZH = 1;
  }
  else {
    WH = 0;
    qH = 0;
    gH = 0;
    HH = 0;
    pH = 0;
    ZH = 0;
  }
  if (graviton()) {
    WG = 1 * 3./4.;
    qG = 3 * 2./3.;
    gG = 1 * 3./4.;
    pG = 1 * 3./4.;
    ZG = 1 * 3./4.;
    GG = 1;
  }
  else {
    WG = 0;
    qG = 0;
    gG = 0;
    pG = 0;
    ZG = 0;
    GG = 0;
  }

  double uu = 6;
  double gW = 1;
  double Wp = 1;
  double WZ = 1;
  double qg = 3;
  double qW = 3;
  double qp = 3;
  double qZ = 3;
  double gp = 1;
  double WW = 1;
  double ZZ = 1;

  double ud = 9;
  double gg = 1;
  double gZ = 1;
  double pp = 1;
  double pZ = 1;

  double tot;
  double dof[18];

  // uu (Q=4/3,I=0) 3 X 3 = 6 + 3-bar
  dof[0]  = 6./9. * uu ;
  dof[0] += 3./9. * uu;
  dof[1]  = 3./9. * ql;
  tot = 0.0;
  for (int i=0;i<2;i++) {
    tot += dof[i];
  }
  br40[0] = dof[0] / tot;
  for (int i=1;i<2;i++) {
    br40[i] = br40[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=4/3,I=0)\n");
    printf("uu %5.1f %f\n",100.0*dof[0]/tot,br40[0]);
    printf("dl %5.1f %f\n",100.0*dof[1]/tot,br40[1]);
  }

  // ud-bar (Q=1,I=0) 3 x 3-bar = 8 + 1
  dof[0]  = 8./9. * ud;
  dof[0] += 1./9. * ud;
  dof[1]  = 8./9. * gW;
  dof[2]  = 1./9. * nl;
  dof[3]  = 1./9. * Wp;
  dof[4]  = 1./9. * WZ;
  dof[5]  = 1./9. * WH;
  dof[6]  = 1./9. * WG;
  tot = 0.0;
  for (int i=0;i<7;i++) {
    tot += dof[i];
  }
  br30[0] = dof[0] / tot;
  for (int i=1;i<7;i++) {
    br30[i] = br30[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=1,I=0)\n");
    printf("ud %5.1f %f\n",100.0*dof[0]/tot,br30[0]);
    printf("gW %5.1f %f\n",100.0*dof[1]/tot,br30[1]);
    printf("nl %5.1f %f\n",100.0*dof[2]/tot,br30[2]);
    printf("Wp %5.1f %f\n",100.0*dof[3]/tot,br30[3]);
    printf("WZ %5.1f %f\n",100.0*dof[4]/tot,br30[4]);
    printf("WH %5.1f %f\n",100.0*dof[5]/tot,br30[5]);
    printf("WG %5.1f %f\n",100.0*dof[6]/tot,br30[6]);
  }

  // ug (Q=2/3,I=1)  3 X 8 = 3 + 6-bar + 15
  dof[0]  = 15./24. * qg;
  dof[0] +=  6./24. * qg;
  dof[0] +=  3./24. * qg;
  dof[1]  =  3./24. * qW;
  dof[2]  =  3./24. * qp;
  dof[3]  =  3./24. * qZ;
  dof[4]  =  3./24. * qH;
  dof[5]  =  3./24. * qG;
  tot = 0.0;
  for (int i=0;i<6;i++) {
     tot += dof[i];
  }
  br21[0] = dof[0] / tot;
  for (int i=1;i<6;i++) {
    br21[i] = br21[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=2/3,I=1)\n");
    printf("ug %5.2f %f\n",100.0*dof[0]/tot,br21[0]);
    printf("dW %5.2f %f\n",100.0*dof[1]/tot,br21[1]);
    printf("up %5.2f %f\n",100.0*dof[2]/tot,br21[2]);
    printf("uZ %5.2f %f\n",100.0*dof[3]/tot,br21[3]);
    printf("uH %5.2f %f\n",100.0*dof[4]/tot,br21[4]);
    printf("uG %5.2f %f\n",100.0*dof[5]/tot,br21[5]);
  }

  // d-bar d-bar (Q=2/3,I=0) 3-bar X 3-bar = 6-bar + 3
  dof[0]  = 6./9. * uu;
  dof[0] += 3./9. * uu;
  dof[1]  = 3./9. * qn;
  dof[2]  = 3./9. * ql;
  tot = 0.0;
  for (int i=0;i<3;i++) {
    tot += dof[i];
  }
  br20[0] = dof[0] / tot;
  for (int i=1;i<3;i++) {
    br20[i] = br20[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=2/3,I=0)\n");
    printf("dd %5.2f %f\n",100.0*dof[0]/tot,br20[0]);
    printf("un %5.2f %f\n",100.0*dof[1]/tot,br20[1]);
    printf("dl %5.2f %f\n",100.0*dof[2]/tot,br20[2]);
  }

  // d-bar g (Q=1/3,I=1) 3-bar X 8 = 3-bar + 6 + 15-bar
  dof[0]  = 15./24. * qg;
  dof[0] +=  6./24. * qg;
  dof[0] +=  3./24. * qg;
  dof[1]  =  3./24. * qW;
  dof[2]  =  3./24. * qp;
  dof[3]  =  3./24. * qZ;
  dof[4]  =  3./24. * qH;
  dof[5]  =  3./24. * qG;
  tot = 0.0;
  for (int i=0;i<6;i++) {
     tot += dof[i];
  }
  br11[0] = dof[0] / tot;
  for (int i=1;i<6;i++) {
    br11[i] = br11[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=1/3,I=1)\n");
    printf("dg %5.2f %f\n",100.0*dof[0]/tot,br11[0]);
    printf("uW %5.2f %f\n",100.0*dof[1]/tot,br11[1]);
    printf("dp %5.2f %f\n",100.0*dof[2]/tot,br11[2]);
    printf("dZ %5.2f %f\n",100.0*dof[3]/tot,br11[3]);
    printf("dH %5.2f %f\n",100.0*dof[4]/tot,br11[4]);
    printf("dG %5.2f %f\n",100.0*dof[5]/tot,br11[5]);
  }

  // ud (Q=1/3,I=0) 3 x 3 = 6 + 3-bar
  dof[0]  = 6./9. * ud;
  dof[0] += 3./9. * ud;
  dof[1]  = 3./9. * qn;
  dof[2]  = 3./9. * ql;
  tot = 0.0;
  for (int i=0;i<3;i++) {
     tot += dof[i];
  }
  br10[0] = dof[0] / tot;
  for (int i=1;i<3;i++) {
    br10[i] = br10[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=1/3,I=0)\n");
    printf("ud %5.2f %f\n",100.0*dof[0]/tot,br10[0]);
    printf("dn %5.2f %f\n",100.0*dof[1]/tot,br10[1]);
    printf("ul %5.2f %f\n",100.0*dof[2]/tot,br10[2]);
  }

  // q q-bar (Q=0,I=0) 3 x 3-bar = 8 + 1
  dof[0]   = 8./9. * ud;
  dof[0]  += 8./9. * ud;
  dof[0]  += 1./9. * ud;
  dof[0]  += 1./9. * ud;
  dof[1]   = 8./9. * gg;
  dof[2]   = 8./9. * gp;
  dof[3]   = 8./9. * gZ;
  dof[4]   = 1./9. * ll;
  dof[5]   = 1./9. * nn;
  dof[6]   = 1./9. * WW;
  dof[7]   = 1./9. * pp;
  dof[8]   = 1./9. * ZZ;
  dof[9]   = 1./9. * pZ;
  dof[10]  = 8./9. * gH;
  dof[11]  = 1./9. * pH;
  dof[12]  = 1./9. * ZH;
  dof[13]  = 1./9. * HH;
  dof[14]  = 8./9. * gG;
  dof[15]  = 1./9. * pG;
  dof[16]  = 1./9. * ZG;
  dof[17]  = 1./9. * GG;
  tot = 0.0;
  for (int i=0;i<18;i++) {
     tot += dof[i];
  }
  br00[0] = dof[0] / tot;
  for (int i=1;i<18;i++) {
    br00[i] = br00[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=0,I=0)\n");
    printf("uu %5.2f %f\n",100.0*dof[0] /tot,br00[0]/2.0);
    printf("dd %5.2f %f\n",100.0*dof[0] /tot,br00[0]/2.0);
    printf("gg %5.2f %f\n",100.0*dof[1] /tot,br00[1] );
    printf("gp %5.2f %f\n",100.0*dof[2] /tot,br00[2] );
    printf("gZ %5.2f %f\n",100.0*dof[3] /tot,br00[3] );
    printf("ll %5.2f %f\n",100.0*dof[4] /tot,br00[4] );
    printf("nn %5.2f %f\n",100.0*dof[5] /tot,br00[5] );
    printf("WW %5.2f %f\n",100.0*dof[6] /tot,br00[6] );
    printf("pp %5.2f %f\n",100.0*dof[7] /tot,br00[7] );
    printf("ZZ %5.2f %f\n",100.0*dof[8] /tot,br00[8] );
    printf("pZ %5.2f %f\n",100.0*dof[9] /tot,br00[9] );
    printf("gH %5.2f %f\n",100.0*dof[10]/tot,br00[10]);
    printf("pH %5.2f %f\n",100.0*dof[11]/tot,br00[11]);
    printf("ZH %5.2f %f\n",100.0*dof[12]/tot,br00[12]);
    printf("HH %5.2f %f\n",100.0*dof[13]/tot,br00[13]);
    printf("gG %5.2f %f\n",100.0*dof[14]/tot,br00[14]);
    printf("pG %5.2f %f\n",100.0*dof[15]/tot,br00[15]);
    printf("ZG %5.2f %f\n",100.0*dof[16]/tot,br00[16]);
    printf("GG %5.2f %f\n",100.0*dof[17]/tot,br00[17]);
  }

  ud = 9 * 5./12.;
  gg = 1;
  gZ = 1;
  pp = 1;
  pZ = 1;
  if (sm()) {
    ll = 3 * 5./12.;
    if (chiral()) {
      nn = 3;
    }
    else {
      nn = 3 * 1./4.;
    }
  }
  else {
    ll =  9 * 5./12.;
    if (majorana()) {
      if (chiral()) {
        nn = 6 * 5./12.;
      }
      else {
        nn = 6 * 1./6.;
      }
    }
    else {
      if (chiral()) {
        nn = 21 * 5./12.;
      }
      else {
        nn = 21 * 1./4.;
      }
    }

  }
  if (higgs()) {
    ZH = 1 * 1./4.;
  }
  else {
    ZH = 0;
  }
  if (graviton()) {
    GG = 1;
  }
  else {
    GG = 0;
  }
  if (graviton()) {
    HG = 1 * 7./12.;
  }
  else {
    HG = 0;
  }
  // gg (Q=0,I=2) 8 x 8 = 1 + 8 + 8 + 10 + 10-bar + 27
  dof[0]   = 16./64. * ud;
  dof[0]  += 16./64. * ud;
  dof[0]  +=  1./64. * ud;
  dof[0]  +=  1./64. * ud;
  dof[1]   = 47./64. * gg;
  dof[1]  += 16./64. * gg;
  dof[1]  +=  1./64. * gg;
  dof[2]   = 16./64. * gp;
  dof[3]   = 16./64. * gZ;
  dof[4]   =  1./64. * ll;
  dof[5]   =  1./64. * nn;
  dof[6]   =  1./64. * WW;
  dof[7]   =  1./64. * pp;
  dof[8]   =  1./64. * ZZ;
  dof[9]   =  1./64. * pZ;
  dof[10]  =  1./64. * ZH;
  dof[11]  =  1./64. * HH;
  dof[12]  =  1./64. * HG;
  dof[13]  =  1./64. * GG;
  tot = 0.0;
  for (int i=0;i<14;i++) {
     tot += dof[i];
  }
  br02[0] = dof[0] / tot;
  for (int i=1;i<14;i++) {
    br02[i] = br02[i-1] + dof[i] / tot;
  }
  if (output) {
    printf("\n(Q=0,I=2)\n");
    printf("uu %5.2f %f\n",100.0*dof[0] /tot,br02[0]/2.0);
    printf("dd %5.2f %f\n",100.0*dof[0] /tot,br02[0]/2.0);
    printf("gg %5.2f %f\n",100.0*dof[1] /tot,br02[1] );
    printf("gp %5.2f %f\n",100.0*dof[2] /tot,br02[2] );
    printf("gZ %5.2f %f\n",100.0*dof[3] /tot,br02[3] );
    printf("ll %5.2f %f\n",100.0*dof[4] /tot,br02[4] );
    printf("nn %5.2f %f\n",100.0*dof[5] /tot,br02[5] );
    printf("WW %5.2f %f\n",100.0*dof[6] /tot,br02[6] );
    printf("pp %5.2f %f\n",100.0*dof[7] /tot,br02[7] );
    printf("ZZ %5.2f %f\n",100.0*dof[8] /tot,br02[8] );
    printf("pZ %5.2f %f\n",100.0*dof[9] /tot,br02[9] );
    printf("ZH %5.2f %f\n",100.0*dof[10]/tot,br02[10]);
    printf("HH %5.2f %f\n",100.0*dof[11]/tot,br02[11]);
    printf("HG %5.2f %f\n",100.0*dof[12]/tot,br02[12]);
    printf("GG %5.2f %f\n",100.0*dof[13]/tot,br02[13]);
  }

  return;
}
/*****************************************************************************/
/* Calculate parton-level cross section for random x_min.                    */
/* Return: xmin   = random shat/s (minimum x value).                         */
/*         Q      = QCD scale (GeV) for parton distribution functions.       */
/*         sighat = parton cross section (pb) at random xmin.                */
/*****************************************************************************/
void QuantumBlackHole::xpart(double& xmin,double& Q,double& sighat)
{
  const double GEV2PB = 389379.0e3;
  const double PI     = acos(-1.0);
  const double pois[] = {0.74,0.61,0.54,0.51,0.50,0.51,0.54};
  const double yrf[]  = {1.54,2.15,2.52,2.77,2.95,3.09,3.20};

  static bool first = true;
  static double rhfact, bhpow, plpow, rpow, a0, a1;
  if (first) {

    // If Randall-Sundrum, set dimensions and Planck scale definition.
    if (RS1()) {
      setTotdim(5);
      setPlanckdef(1);
    }

    // Use Giddings-Thomas Planck scale defintiion in calculations.
    double intmpl; 
    // Convert mplanck to Giddings-Thomas if Randall-Sundrum definitions.
    if (planckdef() == 1) {
      double rs12gt = pow(4.0*PI,1.0/3.0);
      intmpl = rs12gt * mplanck();
    }
    // Convert mplanck to Giddings-Thomas if Dimopoulos-Landsberg..
    else if (planckdef() == 2) {
      double dl2gt = pow( pow(2.0,totdim()-6)*pow(PI,totdim()-5),
                          1.0/(double)(totdim()-2) );
      intmpl = dl2gt * mplanck();
    }
    // Convert mplanck to Giddings-Thomas if PDG definitions.
    else if (planckdef() == 3) {
      double pdg2gt = pow(2,1.0/(double)(totdim()-2));
      intmpl = pdg2gt * mplanck();
    }
    // Leave in Giddings-Thomas definition.
    else {
      intmpl = mplanck();
    }

    // Calcualte other parameters.
    plpow = 2.0 / (double)(totdim()-3);
    // Attempt to increase efficiency and accuracy by power-law transformation.
    int q = qstate();
    int i = istate();
    //bhpow = plpow - 7.0;
    if (q == 4) {
      bhpow = plpow - 2.00665 + 1.0;
    }
    else if (q == -4) {
      bhpow = plpow - 4.24020 + 1.0;
    }
    else if (q == 3) {
      bhpow = plpow - 3.06313 + 1.0;
    }
    else if (q == -3) {
      bhpow = plpow - 3.45916 + 1.0;
    }
    else if (q == 2 && i == 0) {
      bhpow = plpow - 4.09413 + 1.0;
    }
    else if (q == 2 && i == 1) {
      bhpow = plpow - 3.38092 + 1.0;
    }
    else if (q == -2 && i == 0) {
      bhpow = plpow - 2.55985 + 1.0;
    }
    else if (q == -2 && i == 1) {
      bhpow = plpow - 4.28119 + 1.0;
    }
    else if (q == 1 && i == 0) {
      bhpow = plpow - 2.28890 + 1.0;
    }
    else if (q == 1 && i == 1) {
      bhpow = plpow - 4.21862 + 1.0;
    }
    else if (q == -1 && i == 0) {
      bhpow = plpow - 4.16581 + 1.0;
    }
    else if (q == -1 && i == 1) {
      bhpow = plpow - 3.59395 + 1.0;
    }
    else if (q == 0 && i == 0) {
      bhpow = plpow - 3.18361 + 1.0;
    }
    else if (q == 0 && i == 2) {
      bhpow = plpow - 4.34926 + 1.0;
    }
    else {
      //bhpow = plpow - 8.0 + 1.0;
      bhpow = plpow - 2.5978 + 1.0;
    }
    rpow = 1.0 /  bhpow;

    // Area of unit (D-2)-sphere.
    double omegad; 
    if (totdim() == 11) {
      omegad = pow(PI,5) / 12.0;
    }
    else if (totdim() == 10) {
      omegad = 32.0 * pow(PI,4) / 105.0;
    }
    else if (totdim() == 9) {
      omegad = pow(PI,4) / 3.0;
    }
    else if (totdim() == 8) {
      omegad = 16.0 * pow(PI,3) / 15.0;
    }
    else if (totdim() == 7) {
      omegad = pow(PI,3);
    }
    else if (totdim() == 6) {
      omegad = 8.0 * pow(PI,2) / 3.0;
    }
    else if (totdim() == 5) {
      omegad = 2.0 * PI * PI;
    }
    else {
      printf("QBH::xpart: Invalid value for totdim %d\n",totdim());
      printf("            Must be between 5 and 11\n");
      exit(EXIT_FAILURE);
    }
    double mpfact = pow(intmpl,totdim()-2);
    rhfact = 4.0 * pow(2.0*PI,totdim()-4)
           / ((totdim() - 2) * omegad * mpfact);
    first = false;
  }

  // Calculate inelasisity.
  double z = 0.5;
  if (trap()) {
    z = Random::ran();
    int index = -9;
    for (int i=1;i<fsize[totdim()-4];i++) {
      if (za[i] > z && za[i-1] <= z) {
        index = i;
      }
    } 
    if (index < 0) {
      printf("QBH::xpart: Bad index %d %f\n",index,z);
      exit(1);   
    }
    double m = (ya[index] - ya[index-1]) / (za[index] - za[index-1]);
    double b = (ya[index-1]*za[index] - ya[index]*za[index-1]) 
             / (za[index] - za[index-1]);
    double y = m*z + b;
    if (y < 0.0 || y > 1.0) {
      printf("QBH;;xpart: Bad inelasistity y %f\n",y);
      exit(1);
    }
    if (y < minmass()/maxmass()) {
      a0 = pow(maxmass(),bhpow);
      a1 = 0.0;
    }
    else {
      a0 = pow(minmass()/y,bhpow);
      a1 = pow(maxmass(),bhpow) - a0;
    }
  }
  else {
    a0 = pow(minmass(),bhpow);
    a1 = pow(maxmass(),bhpow) - a0;
  }
  
  // Select a random mass for the produced black hole.
  double sqrt_s    = ecm();
  double sqrt_shat = pow(a0 + a1*Random::ran(),rpow);
  double shat      = sqrt_shat * sqrt_shat;
  double emj       = pow(sqrt_shat,1.0-bhpow) / bhpow * a1;

  // Calculate r_h^2.
  double rho2 = pow(rhfact*(sqrt_shat),plpow);

  // Set minimum momentum fraction.
  xmin = shat / (sqrt_s*sqrt_s);

  // Choose QCD momentum scale for parton distribution functions.
  if (qscale()) {
    Q = 1.0 / sqrt(rho2);
  }
  else {
    Q = sqrt_shat;
  }

  // Calculate parton cross section.
  if (z < 0.0 || z > 1.0) {
    printf("bad z %f\n",z);
    exit(1);
  }
  sighat = -GEV2PB * (2.0*z) * PI * rho2 * emj * log(xmin) * 2.0 / (sqrt_shat);

  // Include Poisson probabilites in cross section.
  if (poisson()) sighat *= pois[totdim()-5];

  // Include Yoshino-Rychkov factor.
  if (yrform()) sighat *= yrf[totdim()-5];

  return;
}
/*****************************************************************************/
/* Get the parton distribution functions via PYTHIA PDF class.               */
/* Currently only works for proton-proton collisions.                        */
/* Input:  xx[2]  = x values for both partons.                               */
/*         Q      = QCD scale for parton distribution functions (GeV).       */
/*         pdfPtr = point to PDF class (warning: global variable).           */
/* Return: disf[2][12] = x*pdf(x,Q2) values (warning: class data member) .   */
/*****************************************************************************/
void QuantumBlackHole::pdf(const double* xx,const double Q)
{
  const double q2 = Q*Q;

  for (int i=0;i<6;i++) {
    disf[0][i]   = pdfPtr->xf( i+1,xx[0],q2);
    disf[0][i+6] = pdfPtr->xf(-i-1,xx[0],q2);
  }
  disf[0][12] = pdfPtr->xf(21,xx[0],q2);

  for (int i=0;i<6;i++) {
    disf[1][i]   = pdfPtr->xf( i+1,xx[1],q2);
    disf[1][i+6] = pdfPtr->xf(-i-1,xx[1],q2);
  }
  disf[1][12] = pdfPtr->xf(21,xx[1],q2);

  return;
}
/*****************************************************************************/
/* Sample the cross section to determine the beam partons and return weight. */
/* Input:  sighat = parton-parton cross section (pb).                        */
/*         dxsec  = differental cross section (pb).                          */
/* Return: dxsec  = proton-proton cross section (pb) used as event weight.   */
/*         idn    = PDG ID codes of two partons causing reaction.            */
/*****************************************************************************/
bool QuantumBlackHole::sample(const double sighat,double& dxsec,int* idn)
{
  int ip=-1, iq=-1;

  int q  = qstate();
  int in = istate();
  double rcs = dxsec * Random::ran();
  double sxsec = 0.0;

  // charge 0, qq or gg.
  if (q == 0) {
    if (in == 0) {
      for (int i=1;i<6;i+=2) {
        for (int j=7;j<12;j+=2) {
	  sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done; 
          }
        }
      }
      for (int i=1;i<6;i+=2) {
        for (int j=7;j<12;j+=2) {
	  sxsec += sighat * disf[0][j] * disf[1][i];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done; 
          }
        }
      }
      for (int i=0;i<5;i+=2) {
        for (int j=6;j<11;j+=2) {
	  sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
      for (int i=0;i<5;i+=2) {
        for (int j=6;j<11;j+=2) {
	  sxsec += sighat * disf[0][j] * disf[1][i];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
    }
    if (in == 2) {
      sxsec += sighat * disf[0][12] * disf[1][12];
      if (sxsec > rcs) {
        ip = 12;
        iq = 12;
        goto done;
      }
    }
  }
  // charge 1/3, qq or qg.
  else if (q == 1) {
    if (in == 0) {
      for (int i=1;i<6;i+=2) {
        for (int j=0;j<5;j+=2) {
          sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
      for (int i=1;i<6;i+=2) {
        for (int j=0;j<5;j+=2) {
          sxsec += sighat * disf[0][j] * disf[1][i];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
    }
    if (in == 1) {
      for (int i=6;i<11;i+=2) {
        sxsec += sighat * disf[0][12] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
      for (int i=6;i<11;i+=2) {
        sxsec += sighat * disf[0][i] * disf[1][12];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
    }
  }
  // charge -1/3, qq or qg.
  else if (q == -1) {
    if (in == 0) {
      for (int i=7;i<12;i+=2) {
        for (int j=6;j<11;j+=2) {
          sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
      for (int i=7;i<12;i+=2) {
        for (int j=6;j<11;j+=2) {
          sxsec += sighat * disf[0][j] * disf[1][i];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
    }
    if (in == 1) {
      for (int i=0;i<5;i+=2) {
        sxsec += sighat * disf[0][12] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
      for (int i=0;i<5;i+=2) {
        sxsec += sighat * disf[0][i] * disf[1][12];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
    }
  }
  // charge 2/3, qq or qg.
  else if (q == 2) {
    if (in == 0) {
      for (int i=6;i<11;i+=2) {
        for (int j=6;j<11;j+=2) {
          sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
    }
    if (in == 1) {
      for (int i=1;i<6;i+=2) {
        sxsec += sighat * disf[0][12] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
      for (int i=1;i<6;i+=2) {
        sxsec += sighat * disf[0][i] * disf[1][12];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
    }
  }
  // charge -2/3, qq or qg.
  else if (q == -2) {
    if (in == 0) {
      for (int i=0;i<5;i+=2) {
        for (int j=0;j<5;j+=2) {
          sxsec += sighat * disf[0][i] * disf[1][j];
          if (sxsec > rcs) {
            ip = i;
            iq = j;
            goto done;
          }
        }
      }
    }
    if (in == 1) {
      for (int i=7;i<12;i+=2) {
        sxsec += sighat * disf[0][12] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
      for (int i=7;i<12;i+=2) {
        sxsec += sighat * disf[0][i] * disf[1][12];
        if (sxsec > rcs) {
          ip = i;
          iq = 12;
          goto done;
        }
      }
    }
  }
  // charge 1, qq only.
  else if (q == 3) {
    for (int i=1;i<6;i+=2) {
      for (int j=6;j<11;j+=2) {
        sxsec += sighat * disf[0][i] * disf[1][j];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
    for (int i=1;i<6;i+=2) {
      for (int j=6;j<11;j+=2) {
        sxsec += sighat * disf[0][j] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
  }
  /// charge -1, qq only.
  else if (q == -3) {
    for (int i=0;i<5;i+=2) {
      for (int j=7;j<12;j+=2) {
        sxsec += sighat * disf[0][i] * disf[1][j];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
    for (int i=0;i<5;i+=2) {
      for (int j=7;j<12;j+=2) {
        sxsec += sighat * disf[0][j] * disf[1][i];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
  }
  // charge 4/3, qq only.
  else if (q == 4) {
    for (int i=1;i<6;i+=2) {
      for (int j=1;j<6;j+=2) {
        sxsec += sighat * disf[0][i] * disf[1][j];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
  }
  // charge -4/3, qq only.
  else if (q == -4) {
    for (int i=7;i<12;i+=2) {
      for (int j=7;j<12;j+=2) {
        sxsec += sighat * disf[0][i] * disf[1][j];
        if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
  }
  // All quantum black holes (inclusive).
  else {
    for (int i=0;i<13;i++) {
      for (int j=0;j<13;j++) {
        sxsec += sighat * disf[0][i] * disf[1][j];
	if (sxsec > rcs) {
          ip = i;
          iq = j;
          goto done;
        }
      }
    }
  }

  return false;

done:

  // Convert parton array positions to PDG codes.
  int id[]={++ip,++iq};
  for (int i=0;i<2;i++) {
    // Antiquark.
    if (id[i] > 6 && id[i] <= 12) {
      idn[i] = -id[i] + 6;
    }
    // Gluon.
    else if (id[i] == 13) {
      idn[i] = 21;
    }
    // Quark.
    else {
      idn[i] = id[i];
    }
  }

  return true;
}
/*****************************************************************************/
/* Calculate proton-proton cross section.                                    */
/* Input:  sighat = parton cross section (pb).                               */
/* Return: dxsec  = differental proton-proton cross section (pb).            */
/*****************************************************************************/
void QuantumBlackHole::xsec(const double sighat,double& dxsec)
{
  int q  = qstate();
  int in = istate();
  dxsec = 0.0;

  // charge 0, qq or gg.
  if (q == 0) {
    if (in == 0) {
      for (int i=1;i<6;i+=2) {
        for (int j=7;j<12;j+=2) {
	  dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
      for (int i=1;i<6;i+=2) {
        for (int j=7;j<12;j+=2) {
	  dxsec += sighat * disf[0][j] * disf[1][i];
        }
      }
      for (int i=0;i<5;i+=2) {
        for (int j=6;j<11;j+=2) {
	  dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
      for (int i=0;i<5;i+=2) {
        for (int j=6;j<11;j+=2) {
	  dxsec += sighat * disf[0][j] * disf[1][i];
        }
      }
    }
    if (in == 2) {
      dxsec += sighat * disf[0][12] * disf[1][12];
    }
  }
  // charge 1/3, qq or qg.
  else if (q == 1) {
    if (in == 0) {
      for (int i=1;i<6;i+=2) {
        for (int j=0;j<5;j+=2) {
          dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
      for (int i=1;i<6;i+=2) {
        for (int j=0;j<5;j+=2) {
          dxsec += sighat * disf[0][j] * disf[1][i];
        }
      }
    }
    if (in == 1) {
      for (int i=6;i<11;i+=2) {
        dxsec += sighat * disf[0][12] * disf[1][i];
      }
      for (int i=6;i<11;i+=2) {
        dxsec += sighat * disf[0][i] * disf[1][12];
      }
    }
  }
  // charge -1/3, qq or qg.
  else if (q == -1) {
    if (in == 0) {
      for (int i=7;i<12;i+=2) {
        for (int j=6;j<11;j+=2) {
          dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
      for (int i=7;i<12;i+=2) {
        for (int j=6;j<11;j+=2) {
          dxsec += sighat * disf[0][j] * disf[1][i];
        }
      }
    }
    if (in == 1) {
      for (int i=0;i<5;i+=2) {
        dxsec += sighat * disf[0][12] * disf[1][i];
      }
      for (int i=0;i<5;i+=2) {
        dxsec += sighat * disf[0][i] * disf[1][12];
      }
    }
  }
  // charge 2/3, qq or qg.
  else if (q == 2) {
    if (in == 0) {
      for (int i=6;i<11;i+=2) {
        for (int j=6;j<11;j+=2) {
          dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
    }
    if (in == 1) {
      for (int i=1;i<6;i+=2) {
        dxsec += sighat * disf[0][12] * disf[1][i];
      }
      for (int i=1;i<6;i+=2) {
        dxsec += sighat * disf[0][i] * disf[1][12];
      }
    }
  }
  // charge -2/3, qq or qg.
  else if (q == -2) {
    if (in == 0) {
      for (int i=0;i<5;i+=2) {
        for (int j=0;j<5;j+=2) {
          dxsec += sighat * disf[0][i] * disf[1][j];
        }
      }
    }
    if (in == 1) {
      for (int i=7;i<12;i+=2) {
        dxsec += sighat * disf[0][12] * disf[1][i];
      }
      for (int i=7;i<12;i+=2) {
        dxsec += sighat * disf[0][i] * disf[1][12];
      }
    }
  }
  // charge 1, qq only.
  else if (q == 3) {
    for (int i=1;i<6;i+=2) {
      for (int j=6;j<11;j+=2) {
        dxsec += sighat * disf[0][i] * disf[1][j];
      }
    }
    for (int i=1;i<6;i+=2) {
      for (int j=6;j<11;j+=2) {
        dxsec += sighat * disf[0][j] * disf[1][i];
      }
    }
  }
  /// charge -1, qq only.
  else if (q == -3) {
    for (int i=0;i<5;i+=2) {
      for (int j=7;j<12;j+=2) {
        dxsec += sighat * disf[0][i] * disf[1][j];
      }
    }
    for (int i=0;i<5;i+=2) {
      for (int j=7;j<12;j+=2) {
        dxsec += sighat * disf[0][j] * disf[1][i];
      }
    }
  }
  // charge 4/3, qq only.
  else if (q == 4) {
    for (int i=1;i<6;i+=2) {
      for (int j=1;j<6;j+=2) {
        dxsec += sighat * disf[0][i] * disf[1][j];
      }
    }
  }
  // charge -4/3, qq only.
  else if (q == -4) {
    for (int i=7;i<12;i+=2) {
      for (int j=7;j<12;j+=2) {
        dxsec += sighat * disf[0][i] * disf[1][j];
      }
    }
  }
  // All quantum black holes (inclusive).
  else {
    for (int i=0;i<13;i++) {
      for (int j=0;j<13;j++) {
        dxsec += sighat * disf[0][i] * disf[1][j];
      }
    }
  }

  return;
}
/*****************************************************************************/
/* Quantum black hole differential cross section calculation.                */
/* Return: dxsec = parton cross section (pb) at random x values.             */
/*****************************************************************************/
void QuantumBlackHole::init(double& dxsec)
{
  // Calculate parton-parton cross section at random xmin.
  // Q also calculated.
  double xmin, Q, sighat;
  xpart(xmin,Q,sighat);

  // Generate random parton x values with constraint xmin.
  double xx[2];
  xx[0] = exp(Random::uniform(0.0,log(xmin)));
  xx[1] = xmin / xx[0];

  // Get parton distribution functions.
  // Results stored in class data memeber disf[2][13].
  pdf(xx,Q);

  // Calculate proton-proton cross section.
  xsec(sighat,dxsec);

  return;
}
/*****************************************************************************/
/* Output initialization information.                                        */
/*****************************************************************************/
void QuantumBlackHole::banner(dLHAup* inf)
{
  // Header.
  printf("\n     QBH 3.02 - April 2020\n\n");
  printf("     Simulation of quantum black hole production and decay\n\n");
  printf("     Generator by D.M. Gingrich\n");

  /// ADD or RS1 space.
  if (RS1()) {
      printf("\n     RS1 black hole");
  }
  else {
      printf("\n     ADD black hole");
  }

  // Black hole charge and initial state.
  int q = qstate();
  int s = istate();
  if (q == 4) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (4/3,0)\n");
    }
    else {
      printf("\nQBH::bannder: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 4/3, only state 0 allowed.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == 3) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (1,0)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 1, only state 0 allowed.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == 2) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (2/3,0)\n");
    }
    else if (s == 1) {
      printf("\n     Process: p + p -> QBH (2/3,1)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 2/3, try state 0 or 1.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == 1) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (1/3,0)\n");
    }
    else if (s == 1) {
      printf("\n     Process: p + p -> QBH (1/3,1)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 1/3, try state 0 or 1.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == 0) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (0,0)\n");
    }
    else if (s == 2) {
      printf("\n     Process: p + p -> QBH (0,2)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 0, try state 0 or 2.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == -1) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (-1/3,0)\n");
    }
    else if (s == 1) {
      printf("\n     Process: p + p -> QBH (-1/3,1)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge -1/3, try state 0 or 1.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == -2) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (-2/3,0)\n");
    }
    else if (s == 1) {
      printf("\n     Process: p + p -> QBH (-2/3,1)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge -2/3, try state 0 or 1.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == -3) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (-1,0)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge 1, only state 0 allowed.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (q == -4) {
    if (s == 0) {
      printf("\n     Process: p + p -> QBH (-4/3,0)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for charge -4/3, only state 0 allowed.\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else {
    if (s == 3) {
      printf("\n     Process: p + p -> QBH (all,3)\n");
    }
    else {
      printf("\nQBH::banner: charge %d state %d, state not known.\n",q,s);
      printf("QBH::banner: for all partons, use state 3.\n\n");
      exit(EXIT_FAILURE);
    }
  }

  // Beam parameters.
  printf("\n     Input conditions for this run:\n");
  printf("       Beam 1 (%d) E = %6.3f TeV\n",inf->idBeamA(),
                                              inf->eBeamA()/1000.0);
  printf("       Beam 2 (%d) E = %6.3f TeV\n",inf->idBeamB(),
                                              inf->eBeamB()/1000.0);
  // Basic parameters.
  printf("       Minimum BH mass = %6.3f TeV\n",minmass()/1000.0);
  printf("       Maximum BH mass = %6.3f TeV\n",maxmass()/1000.0);
  printf("       Planck scale    = %6.3f TeV\n",mplanck()/1000.0);
  printf("       Total number of dimensions = %d\n",totdim());

  // Planck scale definition.
  if (planckdef() == 3) {
    printf("       PDG definition of Planck scale\n");
  }
  else {
    printf("       None PDG definition of Planck scale\n");
  }

  // Choice of PDF.
  if (lhaglue() == 10042) {
    printf("       Using PDF set %d CTEQ6L1\n",lhaglue());
  }
  else if (lhaglue() == 10550) {
    printf("       Using PDF set %d CTEQ6.6\n",lhaglue());
  }
  else if (lhaglue() == 10800) {
    printf("       Using PDF set %d CT10\n",lhaglue());
  }
  else if (lhaglue() == 21000) {
    printf("       Using PDF set %d MSTW2008lo\n",lhaglue());
  }
  else if (lhaglue() == 192800) {
    printf("       Using PDF set %d NNPDF21_100\n",lhaglue());
  }
  else if (lhaglue() == 303600) {
    printf("       Using PDF set %d NNPDF31_nnlo_as_0118\n",lhaglue());
  }
  else if (lhaglue() == 0) {
    printf("       Using internal Pythia PDF set %d CTEQ5L\n",lhaglue());
  }
  else {
    printf("       Using PDF set %d\n",lhaglue());
  }

  // Choice of scale.
  if (qscale()) {
    printf("       Using black hole radius as QCD scale\n");
  }
  else {
    printf("       Using black hole mass as QCD scale\n");
  }

  // Choice of Yoshino-Rychkov factors.
  if (trap()) {
    printf("       Yoshinio-Rychkov trapped surface cross section used\n");
  }
  else {
    printf("       Totally inelastic cross section used\n");
  }
  if (yrform()) {
    printf("       Yoshinio-Rychkov factors used\n");
  }
  else {
    printf("       Yoshinio-Rychkov factors not used\n");
  }

  // Poisson two-paticle prbabilities for cross section.
  if (poisson()) {
    printf("       Two-particle decay probability included\n");
  }
  else {
    printf("       Two-particle decay probability not include\n");
  }

  // Global gauge symmetries.
  if (sm()) {
    printf("       Assuming global gauge symmetries conserved\n");
  }
  else {
    printf("       Assuming global gauge symmetries can be violated\n");
  }

  // Choice of outgoing particles.
  if (higgs()) {
    printf("       Include a SM Higgs as particle\n");
  }
  if (graviton()) {
    printf("       Include graviton as particle\n");
  }

  // Handedness of neutrinos.
  if (chiral()) {
    printf("       Assuming neutrinos are left- and right-handed\n");
  }
  else {
    printf("       Assuming neutrinos are only left-handed\n");
  }

  // Dirac or Majorana neutrinos.
  if (majorana()) {
    printf("       Assuming Majorana neutrinos\n");
  }
  else {
    printf("       Assuming Dirac neutrinos\n");
  }

  // Information on the cross section.
  printf("\n");
  printf("     cross section  %e pb\n",inf->xSec(0));
  printf("     error in xsec  %e pb\n",inf->xErr(0));
  printf("     maximum weight %e\n\n",   inf->xMax(0));

  return;
}
/*****************************************************************************/
/* Insert initialization information into an LHEF file.                       */
/*****************************************************************************/
void QuantumBlackHole::trailer(string filename)
{
  // Open Pythia's lhe file.
  ifstream in("LHEFfile.lhe");
  // Open new (QBH) lhe file.
  ofstream QBHfile;
  QBHfile.open(filename.c_str());

  string line;
  // Rewrite first three lines of file.
  (void)getline(in,line);
  QBHfile << line << endl;
  (void)getline(in,line);
  QBHfile << line << endl;
  (void)getline(in,line);
  QBHfile << line << endl;

  // Add new information.
  QBHfile << fixed << setprecision(1);
  QBHfile << "  " << "QBH 3.00 - May 2015 - Parameters:"                                           << endl;
  QBHfile << "  " << "Generator by D.M. Gingrich"                                                   << endl;
  QBHfile << "  "            << minmass()   << " minmass  " << " Minimum BH mass (GeV)"             << endl;
  QBHfile << "  "            << maxmass()   << " maxmass  " << " Maximum BH mass (GeV)"             << endl;
  QBHfile << "  "            << mplanck()   << " mplanck  " << " Planck scale (GeV)"                << endl;
  QBHfile << "  " << setw(6) << RS1()       << " RS1      " << " ADD or RS1 space"                  << endl;
  QBHfile << "  " << setw(6) << qstate()    << " qstate   " << " Black hole charge state"           << endl;
  QBHfile << "  " << setw(6) << istate()    << " istate   " << " Black hole initial state"          << endl;
  QBHfile << "  " << setw(6) << totdim()    << " totdim   " << " Total number of extra dimensions"  << endl;
  QBHfile << "  " << setw(6) << planckdef() << " planckdef" << " Planck scale definition"           << endl;
  QBHfile << "  " << setw(6) << lhaglue()   << " lhaglue  " << " PDF set"                           << endl;
  QBHfile << "  " << setw(6) << qscale()    << " qscale   " << " QCD scale definition"              << endl;
  QBHfile << "  " << setw(6) << trap()      << " trap     " << " Yoshino-Rychkov trapped surface"   << endl;
  QBHfile << "  " << setw(6) << yrform()    << " yrform   " << " Yoshino-Rychkov factors"           << endl;
  QBHfile << "  " << setw(6) << poisson()   << " poisson  " << " Two-particle decay probability"    << endl;
  QBHfile << "  " << setw(6) << sm()        << " sm       " << " Global gauge symmetries conserved" << endl;
  QBHfile << "  " << setw(6) << higgs()     << " higgs    " << " Include a SM Higgs as particle"    << endl;
  QBHfile << "  " << setw(6) << graviton()  << " graviton " << " Include graviton as particle"      << endl;
  QBHfile << "  " << setw(6) << chiral()    << " chiral   " << " Neutrions left- and right-handed"  << endl;
  QBHfile << "  " << setw(6) << majorana()  << " majorana " << " Majorana neutrinos"                << endl;

  // Rewrite reminder of file.
  int i0, i1, i2, i3, i4, i5, i6, i7;
  double d0, d1, d2, d3;
  QBHfile << scientific << setprecision(6);
  while(getline(in,line))
  {
    if (line == "<init>") {
      QBHfile << line << endl;
      in >> i0 >> i1 >> d0 >> d1 >> i2 >> i3 >> i4 >> i5 >> i6 >> i7;
      i6 = 3;
      QBHfile << "  " << i0 << "  " << i1 << "  " << d0 << "  " << d1 << "  " << i2 << "  " << i3 << "  " << i4 << "  " << i5 << "  " << i6 << "  " << i7;
    }
    else if (line == "<event>") {
      QBHfile << line << endl;
      in >> i0 >> i1 >> d0 >> d1 >> d2 >> d3;
      d0 = 1.0;
      d2 = 1.0/137.0;
      d3 = 0.118;
      QBHfile << i0 << " " << i1 << " " << d0 << " " << d1 << " " << d2 << " " << d3;
    }
    else {
      QBHfile << line << endl;
    }
  }

  QBHfile.close();
  return;
}
/*****************************************************************************/
/* Quantum black hole event generation.                                      */
/* Return True: wxsec = parton cross section (pb) at random x values.        */
/*         Q     = QCD scale (GeV) for parton distribution functions.        */
/*         xx    = x values of beam partons.                                 */
/*         idn   = PDG IDs of beam partons.                                  */
/*****************************************************************************/
bool QuantumBlackHole::event(double& wxsec,double& Q,double* xx,int* idn)
{
  // Calculate parton cross section at random xmin.
  // Also return Q.
  double xmin, dxsec;
  xpart(xmin,Q,dxsec);

  // Generate random parton x values with constraint xmin.
  xx[0] = exp(Random::uniform(0.0,log(xmin)));
  xx[1] = xmin / xx[0];

  // Get parton distribution functions.
  // Results stored in private data member disf[2][13].
  pdf(xx,Q);

  // Calculate differential total cross section.
  xsec(dxsec,wxsec);
  if (wxsec == 0.0) {
    //printf("QBH::event: zero cross section\n");
    //idn[0] = -1;
    //idn[1] = -1;
    idn[0] = 2;
    idn[1] = 2;
    //return false;
  }

  // Sample cross section and return true weight and parton IDs.
  bool status = sample(dxsec,wxsec,idn);

  if (!status) {
  //  printf("QBH::event: zero cross section %f\n",wxsec);
  }

  //return status;
  return true;
}
/*****************************************************************************/
/* Determine decay products based on predetermined branching ratios.         */
/* Return:  idn[3:5] = PDG ID codes for decay particles (2 or 3) in event.   */
/*****************************************************************************/
void QuantumBlackHole::decprd(int* idn)
{
  // Calculate quantum black holes state.
  int q  = qstate();
  int in = istate();
  if (in == 3){ 
    q = pdg->particleData.chargeType(idn[0]) 
      + pdg->particleData.chargeType(idn[1]);
    if (idn[0] == 21 && idn[1] == 21) {
      in = 2;
    }
    else if (idn[0] == 21 || idn[1] == 21) {
      in = 1;
    }
    else {
      in = 0;
    }
  }

  int idd3=0, idd4=0, idd5=0;
  double r = Random::ran();
  // charge 0, qq or gg.
  if (q == 0 && in == 0) {
    if (r < br00[0]) {
      if (r < br00[0]/2.0) {
        idd3 =  Random::ranup();
        idd4 = -Random::ranup();
      }
      else {
        idd3 =  Random::randn();
        idd4 = -Random::randn();
      }
    }
    else if (r < br00[1]) {
      idd3 = 21;
      idd4 = 21;
    }
    else if (r < br00[2]) {
      idd3 = 22;
      idd4 = 21;
    }
    else if (r < br00[3]) {
      idd3 = 23;
      idd4 = 21;
    }
    else if (r < br00[4]) {
      if (sm()) {
        int gen = Random::ranlp();
        idd3 = -gen;
        idd4 =  gen;
      }
      else {
        idd3 = -Random::ranlp();
        idd4 =  Random::ranlp();
      }
    }
    else if (r < br00[5]) {
      if (sm()) {
        int gen = Random::rannu();
        if (majorana()) {
          idd3 = gen;
          idd4 = gen;
        }
        else {
          idd3 = -gen;
          idd4 =  gen;
        }
      }
      else {
        if (majorana()) {
          idd3 = Random::rannu();
          idd4 = Random::rannu();
        }
        else {
          idd3 = -Random::rannu();
          idd4 =  Random::rannu();
        }
      }
    }
    else if (r < br00[6]) {
      idd3 =  24;
      idd4 = -24;
    }
    else if (r < br00[7]) {
      idd3 = 22;
      idd4 = 22;
    }
    else if (r < br00[8]) {
      idd3 = 23;
      idd4 = 23;
    }
    else if (r < br00[9]) {
      idd3 = 22;
      idd4 = 23;
    }
    else if (r < br00[10]) {
      idd3 = 25;
      idd4 = 21;
    }
    else if (r < br00[11]) {
      idd3 = 25;
      idd4 = 22;
    }
    else if (r < br00[12]) {
      idd3 = 25;
      idd4 = 23;
    }
    else if (r < br00[13]) {
      idd3 = 25;
      idd4 = 25;
    }
    else if (r < br00[14]) {
      idd3 = 39;
      idd4 = 21;
    }
    else if (r < br00[15]) {
      idd3 = 22;
      idd4 = 39;
    }
    else if (r < br00[16]) {
      idd3 = 23;
      idd4 = 39;
    }
    else {
      idd3 = 39;
      idd4 = 39;
    }
  }
  if (q == 0 && in == 2) {
    if (r < br02[0]) {
      if (r < br02[0]/2.0) {
        idd3 =  Random::ranup();
        idd4 = -Random::ranup();
      }
      else {
        idd3 =  Random::randn();
        idd4 = -Random::randn();
      }
    }
    else if (r < br02[1]) {
      idd3 = 21;
      idd4 = 21;
    }
    else if (r < br02[2]) {
      idd3 = 21;
      idd4 = 22;
    }
    else if (r < br02[3]) {
      idd3 = 21;
      idd4 = 23;
    }
    else if (r < br02[4]) {
      if (sm()) {
        int gen = Random::ranlp();
        idd3 = -gen;
        idd4 =  gen;
      }
      else {
        idd3 = -Random::ranlp();
        idd4 =  Random::ranlp();
      }
    }
    else if (r < br02[5]) {
      if (sm()) {
        int gen = Random::rannu();
        if (majorana()) {
          idd3 = gen;
          idd4 = gen;
        }
        else {
          idd3 = -gen;
          idd4 =  gen;
        }
      }
      else {
        if (majorana()) {
          idd3 = Random::rannu();
          idd4 = Random::rannu();
        }
        else {
          idd3 = -Random::rannu();
          idd4 =  Random::rannu();
        }
      }
    }
    else if (r < br02[6]) {
      idd3 =  24;
      idd4 = -24;
    }
    else if (r < br02[7]) {
      idd3 = 22;
      idd4 = 22;
    }
    else if (r < br02[8]) {
      idd3 = 23;
      idd4 = 23;
    }
    else if (r < br02[9]) {
      idd3 = 22;
      idd4 = 23;
    }
    else if (r < br02[10]) {
      idd3 = 23;
      idd4 = 25;
    }
    else if (r < br02[11]) {
      idd3 = 25;
      idd4 = 25;
    }
    else if (r < br02[12]) {
      idd3 = 25;
      idd4 = 39;
    }
    else {
      idd3 = 39;
      idd4 = 39;
    }
  }
  // charge 1/3, qq or qg.
  else if (q == 1 && in == 0) {
    if (r < br10[0]) {
      idd3 = Random::ranup();
      idd4 = Random::randn();
    }
    else if (r < br10[1]) {
      idd3 =  Random::rannu();
      idd4 = -Random::randn();
      //idd5 =  Random::randn();
      idd5 =  1; // Avoid "Warning in ProcessContainer::constructProcess: unsuitable recoiler found" in Pythia 8.235.
    }
    else {
      idd3 = -Random::ranlp();
      idd4 = -Random::ranup();
      //idd5 =  Random::randn();
      idd5 =  1; // Avoid "Warning in ProcessContainer::constructProcess: unsuitable recoiler found" in Pythia 8.235.
    }
  }
  else if (q == 1 && in == 1) {
    if (r < br11[0]) {
      idd3 = -Random::randn();
      idd4 =  21;
    }
    else if (r < br11[1]) {
      idd3 =  24;
      idd4 = -Random::ranup();
    }
    else if (r < br11[2]) {
      idd3 =  22;
      idd4 = -Random::randn();
    }
    else if (r < br11[3]) {
      idd3 =  23;
      idd4 = -Random::randn();
    }
    else if (r < br11[4]) {
      idd3 =  25;
      idd4 = -Random::randn();
    }
    else {
      idd3 =  39;
      idd4 = -Random::randn();
    }
  }
  // charge -1/3, qq or qg.
  else if (q == -1 && in == 0) {
    if (r < br10[0]) {
      idd3 = -Random::ranup();
      idd4 = -Random::randn();
    }
    else if (r < br10[1]) {
      if (majorana()) {
        idd3 =  Random::rannu();
      }
      else {
        idd3 = -Random::rannu();
      }
      idd4 =  Random::randn();
      //idd5 = -Random::randn();
      idd5 =  -1; // Avoid "Warning in ProcessContainer::constructProcess: unsuitable recoiler found" in Pythia 8.235.
    }
    else {
      idd3 =  Random::ranlp();
      idd4 =  Random::ranup();
      //idd5 = -Random::randn();
      idd5 =  -1; // Avoid "Warning in ProcessContainer::constructProcess: unsuitable recoiler found" in Pythia 8.235.
    }
  }
  else if (q == -1 && in == 1) {
    if (r < br11[0]) {
      idd3 = Random::randn();
      idd4 = 21;
    }
    else if (r < br11[1]) {
      idd3 = -24;
      idd4 =  Random::ranup();
    }
    else if (r < br11[2]) {
      idd3 = 22;
      idd4 = Random::randn();
    }
    else if (r < br11[3]) {
      idd3 = 23;
      idd4 = Random::randn();
    }
    else if (r < br11[4]) {
      idd3 = 25;
      idd4 = Random::randn();
    }
    else {
      idd3 = 39;
      idd4 = Random::randn();
    }
  }
  // charge 2/3, qq or qg.
  else if (q == 2 && in == 0) {
    if (r < br20[0]) {
      idd3 = -Random::randn();
      idd4 = -Random::randn();
    }
    else if (r < br20[1]) {
      if (majorana()) {
        idd3 =  Random::rannu();
      }
      else {
        idd3 = -Random::rannu();
      }
      idd4 =  Random::ranup();
      idd5 = -1;
    }
    else {
      idd3 = -Random::ranlp();
      idd4 =  Random::randn();
      idd5 = -1;
    }
  }
  else if (q == 2 && in == 1) {
    if (r < br21[0]) {
      idd3 = Random::ranup();
      idd4 = 21;
    }
    else if (r < br21[1]) {
      idd3 = 24;
      idd4 = Random::randn();
    }
    else if (r < br21[2]) {
      idd3 = 22;
      idd4 = Random::ranup();
    }
    else if (r < br21[3]) {
      idd3 = 23;
      idd4 = Random::ranup();
    }
    else if (r < br21[4]) {
      idd3 = 25;
      idd4 = Random::ranup();
    }
    else {
      idd3 = 39;
      idd4 = Random::ranup();
    }
  }
  // charge -2/3, qq or qg.
  else if (q == -2 && in == 0) {
    if (r < br20[0]) {
      idd3 = Random::randn();
      idd4 = Random::randn();
    }
    else if (r < br20[1]) {
      idd3 =  Random::rannu();
      idd4 = -Random::ranup();
      idd5 =  1;
    }
    else {
      idd3 =  Random::ranlp();
      idd4 = -Random::randn();
      idd5 =  1;
    }
  }
  else if (q == -2 && in == 1) {
    if (r < br21[0]) {
      idd3 = -Random::ranup();
      idd4 =  21;
    }
    else if (r < br21[1]) {
      idd3 = -24;
      idd4 = -Random::randn();
    }
    else if (r < br21[2]) {
      idd3 =  22;
      idd4 = -Random::ranup();
    }
    else if (r < br21[3]) {
      idd3 =  23;
      idd4 = -Random::ranup();
    }
    else if (r < br21[4]) {
      idd3 =  25;
      idd4 = -Random::ranup();
    }
    else {
      idd3 =  39;
      idd4 = -Random::ranup();
    }
  }
  // charge 1, qq only.
  else if (q == 3) {
    if (r < br30[0]) {
      idd3 =  Random::ranup();
      idd4 = -Random::randn();
    }
    else if (r < br30[1]) {
      idd3 = 24;
      idd4 = 21;
    }
    else if (r < br30[2]) {
      if (sm()) {
        int gen = Random::ranlp();
        idd3 =  gen+1;
        idd4 = -gen;
      }
      else {
        idd3 =  Random::rannu();
        idd4 = -Random::ranlp();
      }
    }
    else if (r < br30[3]) {
      idd3 = 24;
      idd4 = 22;
    }
    else if (r < br30[4]) {
      idd3 = 24;
      idd4 = 23;
    }
    else if (r < br30[5]) {
      idd3 = 24;
      idd4 = 25;
    }
    else {
      idd3 = 24;
      idd4 = 39;
    }
  }
  // charge -1, qq only.
  else if (q == -3) {
    if (r < br30[0]) {
      idd3 =  Random::randn();
      idd4 = -Random::ranup();
    }
    else if (r < br30[1]) {
      idd3 = -24;
      idd4 =  21;
    }
    else if (r < br30[2]) {
      if (sm()) {
        int gen = Random::ranlp();
        if (majorana()) {
          idd3 =  (gen+1);
        }
        else {
          idd3 = -(gen+1);
        }
        idd4 =   gen;
      }
      else {
        if (majorana()) {
          idd3 =  Random::rannu();
        }
        else {
          idd3 = -Random::rannu();
        }
        idd4 =  Random::ranlp();
      }
    }
    else if (r < br30[3]) {
      idd3 = -24;
      idd4 =  22;
    }
    else if (r < br30[4]) {
      idd3 = -24;
      idd4 =  23;
    }
    else if (r < br30[5]) {
      idd3 = -24;
      idd4 =  25;
    }
    else {
      idd3 = -24;
      idd4 =  39;
    }
  }
  // charge 4/3, qq only.
  else if (q == 4) {
    if (r < br40[0]) {
      idd3 = Random::ranup();
      idd4 = Random::ranup();
    }
    else {
      idd3 = -Random::ranlp();
      idd4 = -Random::randn();
      idd5 =  1;
    }
  }
  // charge -4/3, qq only.
  else if (q == -4) {
    if (r < br40[0]) {
      idd3 = -Random::ranup();
      idd4 = -Random::ranup();
    }
    else {
      idd3 =  Random::ranlp();
      idd4 =  Random::randn();
      idd5 = -1;
    }
  }

  idn[2] = idd3;
  idn[3] = idd4;
  idn[4] = idd5;

  return;
}
/*****************************************************************************/
/* Set the colour connections for each final-state particle.                 */
/* Input:  idn[5] = PDG IDs codes for all particles (4 or 5) in event.       */
/* Return: icol[5][2] = Two colour connections for each final-state particle.*/
/*****************************************************************************/
void QuantumBlackHole::flow(int* idn,int icol[5][2])
{
  // Clear array, makes everything colourless.
  for (int i=0;i<5;i++) {
    for (int j=0;j<2;j++) {
      icol[i][j] = 0;
    }
  }

  // Set colour type for each particle.
  int coltype[5];
  for (int i=0;i<5;i++) {
     if (idn[i] > 0 && idn[i] <= 6) {
      coltype[i] = 1; 
    }
    else if (idn[i] >= -6 && idn[i] < 0) {
      coltype[i] = -1; 
    }
    else if (idn[i] == 21) {
      coltype[i] = 2; 
    }
    else {
      coltype[i] = 0; 
    }
  }

  // quark-quark.
  if (coltype[0] == 1 && coltype[1] == 1) {
    if (coltype[2] == 0 || coltype[3] == 0) {
      icol[0][0] = 101;
      icol[1][1] = 102;
      icol[3][1] = 102;
      icol[4][0] = 101;
      if (idn[1]%2 == 0) {
        idn[1] = -1;
      }
      else {
        idn[1] = -2;
      }
    }
    else {
      icol[0][0] = 101;
      icol[1][0] = 102;
      icol[2][0] = 101;
      icol[3][0] = 102;
    }
  }
 
  // antiquark-antiquark.
  else if (coltype[0] == -1 && coltype[1] == -1) {
    if (coltype[2] == 0 || coltype[3] == 0) {
      icol[0][1] = 101;
      icol[1][0] = 102;
      icol[3][0] = 102;
      icol[4][1] = 101;
      if (idn[1]%2 == 0) {
        idn[1] = 1;
      }
      else {
        idn[1] = 2;
      }
    }
    else {
      icol[0][1] = 101;
      icol[1][1] = 102;
      icol[2][1] = 101;
      icol[3][1] = 102;
    }
  }

  // quark-antiquark.
  else if (coltype[0] == 1 && coltype[1] == -1) {
    // Colourless final state.
    if (coltype[2] == 0 && coltype[3] == 0) {
      icol[0][0] = 101;
      icol[1][1] = 101;
    }
    // Two gluon final state.
    else if (coltype[2] == 2 && coltype[3] == 2) {
      icol[0][0] = 101;
      icol[1][1] = 101;
      icol[2][0] = 102;
      icol[2][1] = 103;
      icol[3][0] = 103;
      icol[3][1] = 102;
    }
    // gluon plus colourless state.
    else if (coltype[2] == 0 && coltype[3] == 2) {
      icol[0][0] = 101;
      icol[1][1] = 102;
      icol[3][0] = 101;
      icol[3][1] = 102;
    }
    // quark-antiquark final state.
    else {
      icol[0][0] = 101;
      icol[1][1] = 102;
      icol[2][0] = 101;
      icol[3][1] = 102;
    }
  }

  // antiquark-quark.
  else if (coltype[0] == -1 && coltype[1] == 1) {
    // Colourless final state.
    if (coltype[2] == 0 && coltype[3] == 0) {
      icol[0][1] = 101;
      icol[1][0] = 101;
    }
    // Two gluon final state.
    else if (coltype[2] == 2 && coltype[3] == 2) {
      icol[0][1] = 101;
      icol[1][0] = 101;
      icol[2][0] = 102;
      icol[2][1] = 103;
      icol[3][0] = 103;
      icol[3][1] = 102;
    }
    // gluon plus colourless state.
    else if (coltype[2] == 0 && coltype[3] == 2) {
      icol[0][1] = 101;
      icol[1][0] = 102;
      icol[3][0] = 102;
      icol[3][1] = 101;
    }
    // quark-antiquark final state.
    else {
      icol[0][1] = 101;
      icol[1][0] = 102;
      icol[2][0] = 102;
      icol[3][1] = 101;
    }
  }

  // quark-gluon.
  else if (coltype[0] == 1 && coltype[1] == 2) {
    // quark-gluon final state.
    if (coltype[2] == 1 || coltype[3] == 2) {
      icol[0][0] = 101;
      icol[1][0] = 102;
      icol[1][1] = 103;
      icol[2][0] = 101;
      icol[3][0] = 102;
      icol[3][1] = 103;
    }
    // quark plus colour singlet state.
    else {
      icol[0][0] = 101;
      icol[1][0] = 102;
      icol[1][1] = 101;
      icol[3][0] = 102;
    }
  }

  // gluon-quark.
  else if (coltype[0] == 2 && coltype[1] == 1) {
    // quark-gluon final state.
    if (coltype[2] == 1 || coltype[3] == 2) {
      icol[0][0] = 102;
      icol[0][1] = 103;
      icol[1][0] = 101;
      icol[2][0] = 101;
      icol[3][0] = 102;
      icol[3][1] = 103;
    }
    // quark plus colour singlet state.
    else {
      icol[0][0] = 102;
      icol[0][1] = 101;
      icol[1][0] = 101;
      icol[3][0] = 102;
    }
  }

  // antiquark-gluon.
  else if (coltype[0] == -1 && coltype[1] == 2) {
    // antiquark-gluon final state.
    if (coltype[2] == -1 || coltype[3] == 2) {
      icol[0][1] = 101;
      icol[1][0] = 102;
      icol[1][1] = 103;
      icol[2][1] = 101;
      icol[3][0] = 102;
      icol[3][1] = 103;
    }
    // antiquark plus colour singlet state.
    else {
      icol[0][1] = 101;
      icol[1][0] = 101;
      icol[1][1] = 102;
      icol[3][1] = 102;
    }
  }

  // gluon-antiquark-gluon.
  else if (coltype[0] == 2 && coltype[1] == -1) {
    // gluon-antiquark final state.
    if (coltype[2] == -1 || coltype[3] == 2) {
      icol[0][0] = 102;
      icol[0][1] = 103;
      icol[1][1] = 101;
      icol[2][1] = 101;
      icol[3][0] = 102;
      icol[3][1] = 103;
    }
    // antiquark plus colour singlet state.
    else {
      icol[0][0] = 101;
      icol[0][1] = 102;
      icol[1][1] = 101;
      icol[3][1] = 102;
    }
  }

  // gluon-gluon. 
  else if (coltype[0] == 2 && coltype[1] == 2) {
    // 4-gluon vertex.
    if (coltype[2] == 2 && coltype[3] == 2) {
      icol[0][0] = 101;
      icol[0][1] = 102;
      icol[1][0] = 102;
      icol[1][1] = 101;
      icol[2][0] = 103;
      icol[2][1] = 104;
      icol[3][0] = 104;
      icol[3][1] = 103;
    }
    // 3-gluon vertex.
    else if (coltype[2] == 2 && coltype[3] == 0) {
      icol[0][0] = 101;
      icol[0][1] = 102;
      icol[1][0] = 103;
      icol[1][1] = 101;
      icol[2][0] = 103;
      icol[2][1] = 102;
    }
    // colour singlet final state.
    else if (coltype[2] == 0 && coltype[3] == 0) {
      icol[0][0] = 101;
      icol[0][1] = 102;
      icol[1][0] = 102;
      icol[1][1] = 101;
    }
    // quark-antiquark final state.
    else {
      icol[0][0] = 101;
      icol[0][1] = 102;
      icol[1][0] = 102;
      icol[1][1] = 101;
      icol[2][0] = 103;
      icol[3][1] = 103;
    }
  }

  else {
    printf("QBH::flow: Invalid colour state %d %d\n",coltype[0],coltype[1]) ;
    printf("           Invalid idn %d %d\n",idn[0],idn[1]) ;
    exit(EXIT_FAILURE);
  }

 return;
}
/*****************************************************************************/
/* Calculate kinematics of decay.                                            */
/* Input:  idn[5] = PDG IDs codes for all particles (4 or 5) in event.       */
/*         lha = pointer to derived LHAup class (for beam information).      */
/* Return True: p1 = 5-vector of one decay product.                          */
/*         p2 = 5-vector of other decay product.                             */
/*****************************************************************************/
void QuantumBlackHole::deckin(int* idn,dLHAup* lha,double* pcmf,
                              double* p1,double* p2)
{
  // Calculate black hole 2-particle decay kinematics.
  pcmf[0] = 0.0;
  pcmf[1] = 0.0;
  pcmf[2] = lha->x1pdf()*lha->eBeamA() - lha->x2pdf()*lha->eBeamB();
  pcmf[3] = lha->x1pdf()*lha->eBeamA() + lha->x2pdf()*lha->eBeamB();
  pcmf[4] = (pcmf[3] + pcmf[2])*(pcmf[3] - pcmf[2]) 
          - pcmf[0]*pcmf[0] - pcmf[1]*pcmf[1];
  if (pcmf[4] > 0.0) {
    pcmf[4] = sqrt(pcmf[4]);
  }
  else {
    printf("QBH::deckin: Bad black hole mass %f\n",pcmf[4]);
    exit(EXIT_FAILURE);
  }

  // Kinematics of decay particle 1 in black hole rest frame.
  double mb = pdg->particleData.m0(idn[2]); 
  double mc = pdg->particleData.m0(idn[3]); 
  p1[0] = 0.0;
  p1[1] = 0.0;
  p1[2] = sqrt((pcmf[4]*pcmf[4] - (mb + mc)*(mb + mc)) 
        *      (pcmf[4]*pcmf[4] - (mb - mc)*(mb - mc))) / (2.0*pcmf[4]);
  p1[3] = (pcmf[4]*pcmf[4] - mc*mc + mb*mb) / (2.0*pcmf[4]);
  p1[4] = mb;

  // Kinematics of decay particle 2 in black hole rest frame.
  p2[0] = 0.0;
  p2[1] = 0.0;
  p2[2] = -p1[2];
  p2[3] = (pcmf[4]*pcmf[4] - mb*mb + mc*mc) / (2.0*pcmf[4]);
  p2[4] = mc;

  // Determine spin of black hole.
  int sstate;
  // gluon-gluon: spin 0 or 2.
  if (idn[0] == 21 && idn[1] == 21 ) {
    if (Random::logic(1.0/4.0)) {
      sstate = 0;
    }
    else {
      sstate = 4;
    }
  }
  // quark-gluon: spin 1/2 or 3/2.
  else if (idn[0] == 21 || idn[1] == 21) {
    if (Random::logic(1.0/3.0)) {
      sstate = 3;
    }
    else {
      sstate = 1;
    }
  }
  // quark-quark: spin 0 or 1.
  else {
    if (Random::logic(1.0/4.0)) {
      sstate = 2;
    }
    else {
      sstate = 0;
    }
  }

  // 2-body phase-space decay.
  Decay* decay = new Decay;
  if (sstate == 0) {
    decay->twoBody(pcmf,p1,p2,p1[2],true);
  }
  else {
    decay->twoBody(pcmf,p1,p2,p1[2],false);
  }
  delete decay;

  return;
}
/*****************************************************************************/
/* Les Houches interface virtual functions.                                  */
/*****************************************************************************/
/* Les Houches initialization function.                                      */
/*****************************************************************************/
bool dLHAup::setInit(void)
{
  QuantumBlackHole* qbh = new QuantumBlackHole;

  // Beam particles (only protons=2212 allowed) and energies (must be in GeV).
  setBeamA(2212,qbh->ecm()/2.0,0,qbh->lhaglue());
  setBeamB(2212,qbh->ecm()/2.0,0,qbh->lhaglue());

  // Events are weighted on input but PYTHIA produces events with weight 1.
  setStrategy(1);

  // Calculate cross section, error, and maximum weight.
  //const int NSEARCH = 100000000; // 0.007%, 23 min
  //const int NSEARCH =  10000000; // 0.03%,   3 min
  //const int NSEARCH =   1000000; // 0.07%,  14 s
  const int NSEARCH =    100000; // 0.2%,    2 s
  //const int NSEARCH =     10000; // 0.7%,    1 s
  double bhwsum = 0.0;
  double bhwsq  = 0.0;
  double bhwmax = 0.0;
  for (int i=0;i<NSEARCH;i++) {
    double bhweight;
    qbh->init(bhweight);
    bhwsum += bhweight;
    bhwsq  += bhweight * bhweight;
    if (bhweight > bhwmax) bhwmax = bhweight;
  }
  bhwsum /= double(NSEARCH);
  bhwsq   = bhwsq/double(NSEARCH) - bhwsum*bhwsum;
  if (bhwsq < 0.0) {
    bhwsq = 0.0;
  }
  else {
    bhwsq  = sqrt( bhwsq/double(NSEARCH) );
  }

  // Add process, cross section, error, and maximum weight (must be in pb).
  addProcess(1,bhwsum,bhwsq,bhwmax);

  // Output initialization information.
  qbh->banner(this);

  delete qbh;
  return true;
}
/*****************************************************************************/
/* Les Houches event function.                                               */
/*****************************************************************************/
bool dLHAup::setEvent(int iproc)
{
  if (iproc != 1) {
    (void)printf("setEvent: Bad process ID %d\n",iproc);
    exit(EXIT_FAILURE);
  }

  static bool first = true;
  int idn[5], icol[5][2];
  double bhwgt, Q, xx[2], p0[5], p1[5], p2[5];

  QuantumBlackHole* qbh = new QuantumBlackHole;

  // Generate event: weight, Q, x values and ID of beams.
  bool success = qbh->event(bhwgt,Q,xx,idn);
  if (!success) {
    if (first) {
      printf("dLHAup::setEvent: zero cross section\n");
      first = false;
    } 
    delete qbh;
    return false;
  }

  // Set weight of process for this event.
  setProcess(iproc,bhwgt,Q,-1.0,-1.0);
  
  // Determine black hole decay products.
  qbh->decprd(idn);

  // Connect colour flow.
  qbh->flow(idn,icol);

  // These set functions must occur here since flow can change incident parton type.
  // Set PDF information.
  // Used here to make x values available to other functions.
  setPdf(idn[0],idn[1],xx[0],xx[1],Q,0.0,0.0,true);
  // Set the event information.
  setIdX(idn[0],idn[1],xx[0],xx[1]);

  // Calculate kinematics of decay.
  qbh->deckin(idn,this,p0,p1,p2);

  // Add particles to event record.
  addParticle(idn[0],-1,0,0,icol[0][0],icol[0][1],0.0,0.0,
              x1pdf()*eBeamA(),x1pdf()*eBeamA(),0.0,0.0,9);
  addParticle(idn[1],-1,0,0,icol[1][0],icol[1][1],0.0,0.0,
             -x2pdf()*eBeamB(),x2pdf()*eBeamB(),0.0,0.0,9);
  p1[3] = sqrt(p1[4]*p1[4]+(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]));
  addParticle(idn[2],1,1,2,icol[2][0],icol[2][1],p1[0],p1[1],p1[2],p1[3],p1[4],0.0,9);
  p2[3] = sqrt(p2[4]*p2[4]+(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]));
  addParticle(idn[3],1,1,2,icol[3][0],icol[3][1],p2[0],p2[1],p2[2],p2[3],p2[4],0.0,9);
  if (idn[4] != 0) {
    addParticle(idn[4],1,1,2,icol[4][0],icol[4][1],0.0,0.0,0.0,0.0,0.0,0.0,9);
  }

  // Add black hole as documentation particle.
  //pytia8175 addParticle(40,3,1,2,0,0,0.0,0.0,0.0,0.0,p0[4],0.0,9);

  // Calculate quantum black holes state and use it in the name.
  int q = qbh->qstate();
  int i = qbh->istate();
  if (i == 3) {
    q = QuantumBlackHole::pdg->particleData.chargeType(idn[0]) 
      + QuantumBlackHole::pdg->particleData.chargeType(idn[1]);
    if (idn[0] == 21 && idn[1] == 21) {
      i = 2;
    }
    else if (idn[0] == 21 || idn[1] == 21) {
      i = 1;
    }
    else {
      i = 0;
    }
  }
  ostringstream buffq, buffi;
  buffi << i;
  if (abs(q) == 4 || abs(q) == 2 || q == 0) {
    buffq << q/2;
    string title = "QBH(" + buffq.str() + "," + buffi.str() + ")";
    QuantumBlackHole::pdg->particleData.readString("40:name = " + title);
  }
  else {
    buffq << q;
    string title = "QBH(" + buffq.str() + "/2," + buffi.str() + ")";
    QuantumBlackHole::pdg->particleData.readString("40:name = " + title);
  }

  delete qbh;
  return success;
}
/*****************************************************************************/

} // End of namespace QBH.
//EOF
