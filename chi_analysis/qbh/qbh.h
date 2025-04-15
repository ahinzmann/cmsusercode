//TOF
/*****************************************************************************/
/* QBH version 3.02 - Quantum Black Hole Generator                           */
/* D.M. Gingrich 30-April-2022                                               */
/*****************************************************************************/
#ifndef __QBH
#define __QBH

namespace QBH {

/*****************************************************************************/
/* Define derived LHAup class from PYTHIA.                                   */
/* Needed to define two virtual functions.                                   */
/*****************************************************************************/
class dLHAup : public LHAup {
public:
  bool setInit();
  bool setEvent(int);
};

/*****************************************************************************/
/* Define static class Random.                                               */
/*****************************************************************************/
class Random {
public:
  Random(void);
  static Rndm*  random;
  static bool   logic(const double);
  static int    ranup(void);
  static int    randn(void);
  static int    ranlp(void);
  static int    rannu(void);
  static double ran(void);
  static double uniform(const double,const double);
};

/*****************************************************************************/
/* Define class Decay.                                                       */
/*****************************************************************************/
class Decay {
public:
  void twoBody(const double*,double*,double*,const double,const bool);
private:
  double pp[5];
  double r[3][3];
  void   rotate(void);
  void   matrix(const double*,const double,const double);
  void   boost(const double*,double*);
  void   randomAzimuth(const double,double&,double&);
};

/*****************************************************************************/
/* Define main class QuantumBlackHole.                                       */
/*****************************************************************************/
/* Information about black hole, extra dimensions, particle content, etc.    */
/* Should be fixed for entire execution of program.                          */
/* ecm():       Proton-proton CMS energy (GeV).                              */
/* mplanck():   Fundamental Planck scale (GeV).                              */
/* minmass():   Minimum black hole mass (GeV).                               */
/* maxmass():   Maximum black hole mass (GeV).                               */
/* qstate():    3 times electric charge (4,3,2,1,0,-1,-2,-3,-4 allowed).     */
/* istate():    initial state (0= q-q, 1= q-g, 2= g-g, 3= all partons).      */
/* totdim():    Total number of space-time dimensions.                       */
/* planckdef(): Definition of Planck scale.                                  */
/*              = 1:  Randall-Sundram definition.                            */
/*              = 2:  Dimopoulos-Landsberg definition.                       */
/*              = 3:  PDG definition.                                        */
/*              else: Giddings-Thomas definition.                            */
/* lhaglue():   LHA PDF glue code.                                           */
/* qscale():    Definition of QCD scale for parton distributions.            */
/*              = false: quantum black hole mass.                            */
/*              = true:  1/radius.                                           */
/* yrform():    Use Yoshino-Rychkov form factors or not.                     */
/* trap():      Use Yoshino-Rychkov trapped surface cross section or not.    */
/* poisson():   Use two-particle Poission probability.                       */
/* RS1():       Use Randall-Sundrum black hole (type-1).                     */
/* sm():        Conserve global symmetries or not.                           */
/* chiral():    Neutrions are left- and right-handed.                        */
/* majorana():  Neutrions are Majorana or Dirac particles.                   */
/* higgs():     Include a Standard Model Higgs as particle.                  */
/* graviton():  Include gravition as particle.                               */
/*****************************************************************************/
/* disf[2][13] Holds parton distibutions functions.                          */
/* br40[2]  branching ratios for Q=4/3, qq-state.                            */
/* br30[7]  branching ratios for Q=1,   qq-state.                            */
/* br21[6]  branching ratios for Q=2/3, qg-state.                            */
/* br20[3]  branching ratios for Q=2/3, qq-state.                            */
/* br11[6]  brachning ratios for Q=1/3, qg-state.                            */
/* br10[3]  branching ratios for Q=1/3, qq-state.                            */
/* br00[18] branching ratios for Q=0,   qq-state.                            */
/* br02[14] branching ratios for Q=0,   gg-state.                            */
/* fsize[8] number of lines in input data files.                             */
/* za[91]   x-value (z, scaled impact parameter) in input data files.        */
/* ya[91]   y-value (y, inelastic factor) in input data files.               */
/*****************************************************************************/
class QuantumBlackHole
{
public: 
  QuantumBlackHole(void) {};
  QuantumBlackHole(Pythia*,bool);
  static Pythia* pdg;
  // Member functions.
  bool event(double&,double&,double*,int*);
  void init(double&);
  void banner(dLHAup*);
  void trailer(string);
  void decprd(int*);
  void flow(int*,int i[][2]);
  void deckin(int*,dLHAup*,double*,double*,double*);
  void setLHAglue(int);
  void setTrap(bool);
  // Methods.
  inline void   setQscale(bool l)    {saveQscale    = l;};
  inline void   setYRform(bool l)    {saveYRform    = l;};
  inline void   setPoisson(bool l)   {savePoisson   = l;};
  inline void   setRS1(bool l)       {saveRS1       = l;};
  inline void   setSM(bool l)        {saveSM        = l;};
  inline void   setChiral(bool l)    {saveChiral    = l;};
  inline void   setMajorana(bool l)  {saveMajorana  = l;};
  inline void   setHiggs(bool l)     {saveHiggs     = l;};
  inline void   setGraviton(bool l)  {saveGraviton  = l;};
  inline void   setQstate(int i)     {saveQstate    = i;};
  inline void   setIstate(int i)     {saveIstate    = i;};
  inline void   setTotdim(int i)     {saveTotdim    = i;};
  inline void   setPlanckdef(int i)  {savePlanckdef = i;};
  inline void   setEcm(double a)     {saveEcm       = a;};
  inline void   setMplanck(double a) {saveMplanck   = a;};
  inline void   setMinmass(double a) {saveMinmass   = a;};
  inline void   setMaxmass(double a) {saveMaxmass   = a;};
  inline bool   qscale(void)         {return saveQscale;};
  inline bool   yrform(void)         {return saveYRform;};
  inline bool   trap(void)           {return saveTrap;};
  inline bool   poisson(void)        {return savePoisson;};
  inline bool   RS1(void)            {return saveRS1;};
  inline bool   sm(void)             {return saveSM;};
  inline bool   chiral(void)         {return saveChiral;};
  inline bool   majorana(void)       {return saveMajorana;};
  inline bool   higgs(void)          {return saveHiggs;};
  inline bool   graviton(void)       {return saveGraviton;};
  inline int    qstate(void)         {return saveQstate;};
  inline int    istate(void)         {return saveIstate;};
  inline int    totdim(void)         {return saveTotdim;}
  inline int    planckdef(void)      {return savePlanckdef;};
  inline int    lhaglue(void)        {return saveLHAglue;};
  inline double ecm(void)            {return saveEcm;};
  inline double mplanck(void)        {return saveMplanck;};
  inline double minmass(void)        {return saveMinmass;};
  inline double maxmass(void)        {return saveMaxmass;};
private:
  // Static variables.
  static bool   saveQscale, saveYRform,  saveTrap,     savePoisson,   saveRS1;
  static bool   saveSM,     saveChiral,  saveMajorana, saveHiggs,     saveGraviton;
  static int    saveQstate, saveIstate,  saveTotdim,   savePlanckdef, saveLHAglue;
  static double saveEcm,    saveMplanck, saveMinmass,  saveMaxmass;
  static PDF* pdfPtr;
  static double br40[2];
  static double br30[7];
  static double br21[6];
  static double br20[3];
  static double br11[6];
  static double br10[3];
  static double br00[18];
  static double br02[14];
  static double za[91];
  static double ya[91];
  static int fsize[8];
  // Automatic storage.
  double disf[2][13];
  // Member functions.
  void xpart(double&,double&,double&);
  void pdf(const double*,const double);
  void xsec(const double,double&);
  bool sample(const double,double&,int*);
};

} // End of namespace QBH.
#endif
//EOF
