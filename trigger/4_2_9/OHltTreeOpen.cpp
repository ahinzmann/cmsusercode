//////////////////////////////////////////////////////////////////
// OpenHLT definitions
//////////////////////////////////////////////////////////////////

#define OHltTreeOpen_cxx

#include "OHltTree.h"

#include "TVector2.h"
#include "TPRegexp.h"
#include <stdio.h>
#include <sstream>

using namespace std;


typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > PtEtaPhiELorentzVector;

bool isSingleEleTrigger(TString triggerName, vector<double>& thresholdEle, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){
	
  TString patternEle = "(OpenHLT_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?((TrkIso[VLT]+)?))$";

  TPRegexp matchThresholdEle(patternEle);

  if (matchThresholdEle.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternEle).MatchS(triggerName);
      thresholdEle.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      caloId.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isSingleEleWPTrigger(TString triggerName, vector<double>& thresholdEle, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){ 
  
  TString patternEle = "(OpenHLT_Ele([0-9]+)_WP([0-9]+)){1}$"; 

  TPRegexp matchThresholdEle(patternEle); 

  if (matchThresholdEle.MatchB(triggerName)) 
    { 
      TObjArray *subStrL   = TPRegexp(patternEle).MatchS(triggerName); 
      thresholdEle.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof()); 
      caloId.push_back("WP"+((TObjString *)subStrL->At(3))->GetString()); 
      caloIso.push_back("WP"+((TObjString *)subStrL->At(3))->GetString()); 
      trkId.push_back("WP"+((TObjString *)subStrL->At(3))->GetString()); 
      trkIso.push_back("WP"+((TObjString *)subStrL->At(3))->GetString()); 
      delete subStrL; 
 
      return true; 
    } 
  else return false; 
} 


bool isDoubleEleTrigger(TString triggerName, vector<double>& thresholdEle, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){
	
  TString patternEle = "(OpenHLT_DoubleEle([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?((TrkIso[VLT]+)?))$";

  TPRegexp matchThresholdEle(patternEle);

  if (matchThresholdEle.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternEle).MatchS(triggerName);
      thresholdEle.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      caloId .push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isAsymDoubleEleTrigger(TString triggerName, vector<double>& thresholdEle, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){
	
  TString patternEle = "(OpenHLT_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?((TrkIso[VLT]+)?))$";

  TPRegexp matchThresholdEle(patternEle);

  if (matchThresholdEle.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternEle).MatchS(triggerName);
      thresholdEle.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      caloId .push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      thresholdEle.push_back((((TObjString *)subStrL->At(7))->GetString()).Atof());
      caloId .push_back(((TObjString *)subStrL->At(8))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(9))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(10))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(11))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isSingleEleX_HTXTrigger(TString triggerName, vector<double>& thresholds, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){
	
  TString patternEle = "(OpenHLT_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_HT([0-9]+))$";

  TPRegexp matchThresholdEle(patternEle);

  if (matchThresholdEle.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternEle).MatchS(triggerName);
      thresholds.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());//Ele
      caloId.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      thresholds.push_back((((TObjString *)subStrL->At(7))->GetString()).Atof());//HT
      delete subStrL;

      return true;
    }
  else return false;
}


bool isSinglePhotonTrigger(TString triggerName, vector<double>& thresholdPhoton, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso){
	
  TString patternPhoton = "(OpenHLT_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?((Iso[VXLT]+)?))$";

  TPRegexp matchThresholdPhoton(patternPhoton);

  if (matchThresholdPhoton.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternPhoton).MatchS(triggerName);
      thresholdPhoton.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isDoublePhotonTrigger(TString triggerName, vector<double>& thresholdPhoton, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso){
	
  TString patternPhoton = "(OpenHLT_DoublePhoton([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?((Iso[VXLT]+)?))$";

  TPRegexp matchThresholdPhoton(patternPhoton);

  if (matchThresholdPhoton.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternPhoton).MatchS(triggerName);
      thresholdPhoton.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}


bool isAsymDoublePhotonTrigger(TString triggerName, vector<double>& thresholdPhoton, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso){
	
  TString patternPhoton = "(OpenHLT_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?(Iso[VXLT]+)?_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?((Iso[VXLT]+)?))$";

  TPRegexp matchThresholdPhoton(patternPhoton);

  if (matchThresholdPhoton.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(patternPhoton).MatchS(triggerName);
      thresholdPhoton.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      thresholdPhoton.push_back((((TObjString *)subStrL->At(6))->GetString()).Atof());
      r9Id.push_back(((TObjString *)subStrL->At(7))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(8))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(9))->GetString());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isMuXTrigger(TString triggerName, vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_Mu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isIsoMuXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_IsoMu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isIsoMuX_eta2pXTrigger(TString triggerName, vector<double> &thresholds)  
{
 TString pattern = "(OpenHLT_IsoMu([0-9]+)_eta2p([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdEta  = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(2. + thresholdEta/10.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMuX_eta2pXTrigger(TString triggerName, vector<double> &thresholds)  
{
 TString pattern = "(OpenHLT_Mu([0-9]+)_eta2p([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdEta  = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(2. + thresholdEta/10.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMuX_MuXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_Mu([0-9]+)_Mu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL    = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu0 = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdL3Mu1 = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu0);
      thresholds.push_back(thresholdL3Mu1);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMuX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_Mu([0-9]+)_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT   = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMuX_EleX_HTXTrigger(TString triggerName, vector<double>& thresholds, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso){
	
  TString pattern = "(OpenHLT_Mu([0-9]+)_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_HT([0-9])+)$";

  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      thresholds.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());
      thresholds.push_back((((TObjString *)subStrL->At(3))->GetString()).Atof());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(6))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(7))->GetString());
      thresholds.push_back((((TObjString *)subStrL->At(8))->GetString()).Atof());
      delete subStrL;

      return true;
    }
  else return false;
}

bool isJetXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Jet([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isJetXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Jet([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isJetX_NoJetIDTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Jet([0-9]+)_NoJetID){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isDiJetAveXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJetAve([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDiJetAveXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJetAve([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isMeffXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Meff([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMeff = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMeff);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMeffXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Meff([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMeff = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMeff);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTX_CentralJetX_BTagIPTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)_CentralJet([0-9]+)_BTagIP){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdJet = (((TObjString *)subStrL->At(3))->GetString()).Atof(); 
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdJet); 
      return true;
    }
  else
    return false;
}

bool isBTagIP_pfMHTX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  //   TString pattern = "(OpenHLT_BTagIP_PFMHT([0-9]+)_HT([0-9]+)){1}$";
  TString pattern = "(OpenHLT_HT([0-9]+)_CentralJet([0-9]+)_BTagIP_PFMHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof(); 
      double thresholdJet = (((TObjString *)subStrL->At(3))->GetString()).Atof(); 
      double thresholdMET = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdHT); 
      thresholds.push_back(thresholdJet); 
      thresholds.push_back(thresholdMET);
      return true;
    }
  else
    return false;
}

bool ispfMHTX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)_PFMHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof(); 
      double thresholdMET = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdHT); 
      thresholds.push_back(thresholdMET);
      return true;
    }
  else
    return false;
}

bool isMHTXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MHT([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMETXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMET = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool ispfMHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_PFMHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdpfMHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdpfMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDiJetXU_PTXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJet([0-9]+)U_PT([0-9]+)U)$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdPT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdPT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDiJetX_PTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJet([0-9]+)_PT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdPT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdPT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTXU_MHTXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)U_MHT([0-9]+)U)$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTX_MHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)_MHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isPFHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_PFHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPFHT  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdPFHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDiCentralPFJetX_PFMHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiCentralPFJet([0-9]+)_PFMHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPFJet  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdPFMHT  = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdPFJet);
      thresholds.push_back(thresholdPFMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isPFHTX_PFMHTX_OrHTX_OrMHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_PFHT([0-9]+)_PFMHT([0-9]+)_OrHT([0-9]+)_OrMHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPFHT  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdPFMHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdHT    = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      double thresholdMHT   = (((TObjString *)subStrL->At(5))->GetString()).Atof();
      thresholds.push_back(thresholdPFHT);
      thresholds.push_back(thresholdPFMHT);
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}



bool isAlphaTTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_HT([0-9]+)_AlphaT0p([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double HtThreshold = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double BetaThreshold = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(HtThreshold);
      thresholds.push_back(BetaThreshold);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isDoubleJetXU_ForwardBackwardTrigger(
					  TString triggerName,
					  vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoubleJet([0-9]+)U_ForwardBackward){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDoubleJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDoubleJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleJetX_ForwardBackwardTrigger(
					 TString triggerName,
					 vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoubleJet([0-9]+)_ForwardBackward){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDoubleJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDoubleJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isQuadJetXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_QuadJet([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdQuadJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdQuadJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isQuadJetXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_QuadJet([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdQuadJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdQuadJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isQuadJetX_BTagIPTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_QuadJet([0-9]+)_BTagIP){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdQuadJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdQuadJet);
      delete subStrL;
      return true;
    }
  else
    return false;
	
}

bool isQuadJetX_IsoPFTauXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_QuadJet([0-9]+)_IsoPFTau([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdQuadJet  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdIsoPFTau = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdQuadJet);
      thresholds.push_back(thresholdIsoPFTau);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isQuadJetX_IsoPFTauX_PFMHTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_QuadJet([0-9]+)_IsoPFTau([0-9]+)_PFMHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdQuadJet  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdIsoPFTau = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdIsoPFMHT = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdQuadJet);
      thresholds.push_back(thresholdIsoPFTau);
      thresholds.push_back(thresholdIsoPFMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}




bool isR0X_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdR = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.); 
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isR0X_MRX_BTagTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_R0([0-9]+)_MR([0-9]+)_CentralJet([0-9]+)_BTagIP{1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdR = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdJet = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isRMRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_RMR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdRMR = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdRMR);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isPhotonX_R0X_MRXTrigger(TString triggerName, vector<double> &thresholds, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso)
{
	
  TString pattern = "(OpenHLT_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?(Iso[VXLT]+)?_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName)) 
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      double thresholdR = (((TObjString *)subStrL->At(6))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(7))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoublePhotonX_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoublePhoton([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoublePhotonX_CaloIdL_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoublePhoton([0-9]+)_CaloIdL_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isDoublePhotonX_R0X_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoublePhoton([0-9]+)_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoublePhotonX_CaloIdL_R0X_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DoublePhoton([0-9]+)_CaloIdL_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      thresholds.push_back(40.);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isMuX_R0X_MRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Mu([0-9]+)_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isEleX_R0X_MRXTrigger(TString triggerName, vector<double> &thresholds,  vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso)
{
	
  TString pattern = "(OpenHLT_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_R0([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      caloId.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      double thresholdR = (((TObjString *)subStrL->At(7))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(8))->GetString()).Atof();
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdR/100.);
      thresholds.push_back(thresholdMR);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isMuX_MR0XTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Mu([0-9]+)_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdMR);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isEleLX_MR0XTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Ele([0-9]+)_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdMR);
      thresholds.push_back(thresholdEle);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isEleTX_MR0XTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_Ele([0-9]+)_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMR = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdMR);
      thresholds.push_back(thresholdEle);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isR0XTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_R0([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdR = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdR/100.);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMRXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MR([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMR = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMR);
      delete subStrL;
      return true;
    }
  else
    return false;
}



bool isPT12U_XUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_PT12U_([0-9]+)U{1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdPT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isPT12_XTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_PT12_([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPT = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdPT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_DiJetXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_DiJet([0-9]+)U{1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_DiJetXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_DiJet([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_JetXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_DiJet([0-9]+)U{1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_JetXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_Jet([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_DiJetXU_MuXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_DiJet([0-9]+)U_Mu([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdL3Mu = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdL3Mu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isBTagMu_DiJetX_MuXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagMu_DiJet([0-9]+)_Mu([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdL3Mu = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdL3Mu);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isBTagIP_JetXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_BTagIP_Jet([0-9]+){1})$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isSingleIsoTauX_TrkX_METXTrigger(
				      TString triggerName,
				      vector<double> &thresholds)
{
	
  TString pattern =
    "(OpenHLT_SingleIsoTau([0-9]+)_Trk([0-9]+)_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdTau = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdTrk = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdTau);
      thresholds.push_back(thresholdTrk);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isIsoPFTauX_TrkX_METXTrigger(
                                      TString triggerName,
                                      vector<double> &thresholds)
{

  TString pattern =
    "(OpenHLT_IsoPFTau([0-9]+)_Trk([0-9]+)_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdTau = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdTrk = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdTau);
      thresholds.push_back(thresholdTrk);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isCentralJetXU_METXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_CentralJet([0-9]+)U_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdCenJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdCenJet);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isCentralJetX_METXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_CentralJet([0-9]+)_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdCenJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdCenJet);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}



bool isDiJetXU_METXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJet([0-9]+)U_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDiJetX_METXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_DiJet([0-9]+)_MET([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMET = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      thresholds.push_back(thresholdMET);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isMETX_HTXUTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MET([0-9]+)_HTX([0-9]+)U){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMET = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdMET);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMETX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_MET([0-9]+)_HTX([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMET = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdMET);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isExclDiJetXU_HFORTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_ExclDiJet([0-9]+)U_HFOR){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isExclDiJetX_HFORTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_ExclDiJet([0-9]+)_HFOR){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isExclDiJetXU_HFANDTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_ExclDiJet([0-9]+)U_HFAND){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isExclDiJetX_HFANDTrigger(TString triggerName, vector<double> &thresholds)
{
	
  TString pattern = "(OpenHLT_ExclDiJet([0-9]+)_HFAND){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdDiJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdDiJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isL1SingleMuXTrigger(TString triggerName)
{
	
  TString pattern = "(OpenHLT_L1SingleMu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      return true;
    }
  else
    return false;
}

bool isL1DoubleMuXTrigger(TString triggerName)
{
	
  TString pattern = "(OpenHLT_L1DoubleMu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      return true;
    }
  else
    return false;
}

bool isL2SingleMuXTrigger(TString triggerName, vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_L2Mu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isL2DoubleMuXTrigger(TString triggerName, vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_L2DoubleMu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdMu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isL1SingleEGXTrigger(TString triggerName)
{
	
  TString pattern = "(OpenHLT_L1SingleEG([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      return true;
    }
  else
    return false;
}

bool isPhotonX_HTXTrigger(TString triggerName, vector<double> &thresholds, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso)
{
	
  TString pattern = "(OpenHLT_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?(Iso[VXLT]+)?_HT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      double thresholdHT = (((TObjString *)subStrL->At(6))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isPhotonX_MHTXTrigger(TString triggerName, vector<double> &thresholds, vector<TString>& r9Id, vector<TString>& caloId,  vector<TString>& photonIso)
{
	
  TString pattern = "(OpenHLT_Photon([0-9]+)_?(R9Id)?_?(CaloId[VXLT]+)?_?(Iso[VXLT]+)?_MHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      r9Id.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      photonIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      double thresholdMHT = (((TObjString *)subStrL->At(6))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

// bool isPhotonX_CaloIdL_HTXTrigger(
// 				  TString triggerName,
// 				  vector<double> &thresholds)
// {
//   // 2011-03-29 promoted to vX TODO check
//   TString pattern = "(OpenHLT_Photon([0-9]+)_CaloIdL_HT([0-9]+)(_v[0-9]+))$";
//   TPRegexp matchThreshold(pattern);
	
//   if (matchThreshold.MatchB(triggerName))
//     {
//       TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
//       double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
//       double thresholdHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
//       thresholds.push_back(thresholdPhoton);
//       thresholds.push_back(thresholdHT);
//       delete subStrL;
//       return true;
//     }
//   else
//     return false;
// }

bool isPhotonX_CaloIdL_MHTXTrigger(
				   TString triggerName,
				   vector<double> &thresholds)
{
  // 2011-03-29 promoted to vX TODO check
  TString pattern = "(OpenHLT_Photon([0-9]+)_CaloIdL_MHT([0-9]+)(_v[0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdPhoton = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdPhoton);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTX_MuX_pfMHTXTrigger( 
			     TString triggerName,
			     vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_HT([0-9]+)_Mu([0-9]+)_PFMHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL    = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT    = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdL3Mu  = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdpfMHT = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdpfMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isEleX_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HTXTrigger(
							 TString triggerName,
							 vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_Ele([0-9]+)_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL    = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle   = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT    = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isEleX_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HTXTrigger(
							  TString triggerName,
							  vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_Ele([0-9]+)_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL    = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle   = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT    = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isHTX_EleX_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_pfMHTXTrigger(
								  TString triggerName,
								  vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_HT([0-9]+)_Ele([0-9]+)_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_PFMHT([0-9]+))$";
  TPRegexp matchThreshold(pattern);
	
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL    = TPRegexp(pattern).MatchS(triggerName);
      double thresholdHT    = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdEle   = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdMHT   = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdMHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isIsoPFTauX_TrkXTrigger( 
                                            TString triggerName, 
                                            vector<double> &thresholds) 
{ 
  TString pattern = "(OpenHLT_IsoPFTau([0-9]+)_Trk([0-9]+))$"; 
  TPRegexp matchThreshold(pattern); 
         
  if (matchThreshold.MatchB(triggerName)) 
    { 
      TObjArray *subStrL     = TPRegexp(pattern).MatchS(triggerName); 
      double thresholdTau0   = (((TObjString *)subStrL->At(2))->GetString()).Atof(); 
      double thresholdTrk    = (((TObjString *)subStrL->At(3))->GetString()).Atof(); 
      thresholds.push_back(thresholdTau0); 
      thresholds.push_back(thresholdTrk); 
      delete subStrL; 
      return true; 
    } 
  else 
    return false; 
} 

bool isDoubleIsoPFTauX_X_TrkX_eta2pXTrigger(
                                            TString triggerName,
                                            vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_DoubleIsoPFTau([0-9]+)_([0-9]+)_Trk([0-9]+)_eta2p([0-9]+))$";
  TPRegexp matchThreshold(pattern);
        
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL     = TPRegexp(pattern).MatchS(triggerName);
      double thresholdTau0   = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdTau1   = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdTrk    = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      double thresholdEta    = (((TObjString *)subStrL->At(5))->GetString()).Atof();
      thresholds.push_back(thresholdTau0);
      thresholds.push_back(thresholdTau1);
      thresholds.push_back(thresholdTrk);
      double eta = 2.+ thresholdEta/10.;
      thresholds.push_back(eta);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleLooseIsoPFTauX_X_TrkX_eta2pXTrigger(
						 TString triggerName,
						 vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_DoubleLooseIsoPFTau([0-9]+)_([0-9]+)_Trk([0-9]+)_eta2p([0-9]+))$";
  TPRegexp matchThreshold(pattern);
        
  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL     = TPRegexp(pattern).MatchS(triggerName);
      double thresholdTau0   = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdTau1   = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdTrk    = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      double thresholdEta    = (((TObjString *)subStrL->At(5))->GetString()).Atof();
      thresholds.push_back(thresholdTau0);
      thresholds.push_back(thresholdTau1);
      thresholds.push_back(thresholdTrk);
      double eta = 2.+ thresholdEta/10.;
      thresholds.push_back(eta);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isFatJetMass(TString triggerName, vector<double> &thresholds)
{

  TString pattern  = "(OpenHLT_FatJetMass([0-9]+)_DR([0-9]+)p([0-9]+)_Deta([0-9]+)p([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMass = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdDR = (((TObjString *)subStrL->At(3))->GetString()).Atof()+(((TObjString *)subStrL->At(4))->GetString()).Atof()/10. ;
      double thresholdDEta = (((TObjString *)subStrL->At(5))->GetString()).Atof()+(((TObjString *)subStrL->At(6))->GetString()).Atof()/10. ;
      thresholds.push_back(thresholdMass);
      thresholds.push_back(thresholdDR);
      thresholds.push_back(thresholdDEta);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isFatJetMassBTag(TString triggerName, vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_FatJetMass([0-9]+)_DR([0-9]+)p([0-9]+)_Deta([0-9]+)p([0-9]+)_CentralJet30_BTagIP){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
      double thresholdMass = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdDR = (((TObjString *)subStrL->At(3))->GetString()).Atof()+(((TObjString *)subStrL->At(4))->GetString()).Atof()/10. ;
      double thresholdDEta = (((TObjString *)subStrL->At(5))->GetString()).Atof()+(((TObjString *)subStrL->At(6))->GetString()).Atof()/10. ;
      //double thresholdJet = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdMass);
      thresholds.push_back(thresholdDR);
      thresholds.push_back(thresholdDEta);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isJetX_BTagIPTrigger(TString triggerName, vector<double> &thresholds)
{
  TString pattern = "(OpenHLT_Jet([0-9]+)_CentralJet30_BTagIP){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL  = TPRegexp(pattern).MatchS(triggerName);
      double thresholdJet = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdJet);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleEleX_CaloIdL_TrkIdVL_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleEle([0-9]+)_CaloIdL_TrkIdVL_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL  = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}
 
bool isDoubleEleX_CaloIdT_TrkIdVL_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleEle([0-9]+)_CaloIdT_TrkIdVL_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL  = TPRegexp(pattern).MatchS(triggerName);
      double thresholdEle = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdEle);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleMuXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleMu([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL  = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleMuX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleMu([0-9]+)_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL   = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdHT   = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleMuX_MassX_HTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleMu([0-9]+)_Mass([0-9]+)_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL      = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu    = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMassCut = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdHT      = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdMassCut);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isDoubleMuX_MassX_HTFJX_PFHTXTrigger(TString triggerName, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleMu([0-9]+)_Mass([0-9]+)_HTFJ([0-9]+)_PFHT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL      = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu    = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      double thresholdMassCut = (((TObjString *)subStrL->At(3))->GetString()).Atof();
      double thresholdHT      = (((TObjString *)subStrL->At(4))->GetString()).Atof();
      double thresholdPFHT      = (((TObjString *)subStrL->At(5))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back(thresholdMassCut);
      thresholds.push_back(thresholdHT);
      thresholds.push_back(thresholdPFHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isDoubleEleX_MassX_HTXTrigger(TString triggerName,  vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso, vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_DoubleEle([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_Mass([0-9]+)_HT([0-9]+)){1}$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL      = TPRegexp(pattern).MatchS(triggerName);
      thresholds.push_back((((TObjString *)subStrL->At(2))->GetString()).Atof());//Ele
      caloId.push_back(((TObjString *)subStrL->At(3))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(4))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(6))->GetString());
      double thresholdMassCut = (((TObjString *)subStrL->At(7))->GetString()).Atof();
      double thresholdHT      = (((TObjString *)subStrL->At(8))->GetString()).Atof();
      thresholds.push_back(thresholdMassCut);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}

bool isMuX_EleX_MassX_HTXTrigger(TString triggerName, vector<TString>& caloId, vector<TString>& caloIso, vector<TString>& trkId, vector<TString>& trkIso,  vector<double> &thresholds)
{
 TString pattern = "(OpenHLT_Mu([0-9]+)_Ele([0-9]+)_?(CaloId[VXLT]+)?_?(CaloIso[VLT]+)?_?(TrkId[VLT]+)?_?(TrkIso[VLT]+)?_Mass([0-9]+)_HT([0-9]+))$";
  TPRegexp matchThreshold(pattern);

  if (matchThreshold.MatchB(triggerName))
    {
      TObjArray *subStrL      = TPRegexp(pattern).MatchS(triggerName);
      double thresholdL3Mu    = (((TObjString *)subStrL->At(2))->GetString()).Atof();
      thresholds.push_back(thresholdL3Mu);
      thresholds.push_back((((TObjString *)subStrL->At(3))->GetString()).Atof());//Ele
      caloId.push_back(((TObjString *)subStrL->At(4))->GetString());
      caloIso.push_back(((TObjString *)subStrL->At(5))->GetString());
      trkId.push_back(((TObjString *)subStrL->At(6))->GetString());
      trkIso.push_back(((TObjString *)subStrL->At(7))->GetString());
      double thresholdMassCut = (((TObjString *)subStrL->At(8))->GetString()).Atof();
      double thresholdHT      = (((TObjString *)subStrL->At(9))->GetString()).Atof();
      thresholds.push_back(thresholdMassCut);
      thresholds.push_back(thresholdHT);
      delete subStrL;
      return true;
    }
  else
    return false;
}


bool isNJetPtTrigger(TString triggerName, vector<double> &thresholds)
{
    TString pattern = "(OpenHLT_([0-9]+)Jet([0-9]+))(_v[0-9]+)?$";
    TPRegexp matchThreshold(pattern);
	
    if (matchThreshold.MatchB(triggerName))
    {
        TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
        double thresholdN  = (((TObjString *)subStrL->At(2))->GetString()).Atof();
        double thresholdPt = (((TObjString *)subStrL->At(3))->GetString()).Atof();
        thresholds.push_back(thresholdN);
        thresholds.push_back(thresholdPt);
        delete subStrL;
        return true;
    }
    else
        return false;
}


bool isNTowCalEt0pTrigger(TString triggerName, vector<double> &thresholds)
{
    TString pattern = "(OpenHLT_([0-9]+)TowCal0p([0-9]+))(_v[0-9]+)?$";
    TPRegexp matchThreshold(pattern);

    if (matchThreshold.MatchB(triggerName))
    {
        TObjArray *subStrL = TPRegexp(pattern).MatchS(triggerName);
        double thresholdN  = ((TObjString *)subStrL->At(2))->GetString().Atof();
        double thresholdEt = 0;
        istringstream iss(string("0.") + ((TObjString *)subStrL->At(3))->GetString().Data());
        iss >> thresholdEt;
        thresholds.push_back(thresholdN);
        thresholds.push_back(thresholdEt);
        delete subStrL;
        return true;
    }
    else
        return false;
}


void OHltTree::CheckOpenHlt(
			    OHltConfig *cfg,
			    OHltMenu *menu,
			    OHltRateCounter *rcounter,
			    int it)
{
  TString triggerName = menu->GetTriggerName(it);
  vector<double> thresholds;
  vector<double> thresholdEle;
  vector<double> thresholdPhoton;
  vector<TString> caloId; 
  vector<TString> caloIso;
  vector<TString> trkId;  
  vector<TString> trkIso; 
  vector<TString> photonIso;
  vector<TString> r9Id;

  //////////////////////////////////////////////////////////////////
  // Check OpenHLT L1 bits for L1 rates
	
  if (triggerName.CompareTo("OpenL1_ZeroBias") == 0)
    {
      if (prescaleResponse(menu, cfg, rcounter, it))
	{
	  triggerBit[it] = true;
	}
    }
  else if (triggerName.CompareTo("OpenL1_EG5_HTT100") == 0)
    {
      if (map_BitOfStandardHLTPath.find(triggerName)->second == 1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenL1_Mu3EG5") == 0)
    {
      if (map_BitOfStandardHLTPath.find(triggerName)->second == 1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenL1_QuadJet8U") == 0)
    {
      if (map_BitOfStandardHLTPath.find(triggerName)->second == 1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
	
  //////////////////////////////////////////////////////////////////
  // Example for pass through triggers
	
  else if (triggerName.CompareTo("OpenHLT_L1Seed1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1Seed2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
	
  //////////////////////////////////////////////////////////////////
  // Check OpenHLT triggers
	

  /*SingleEle*/

  else if (isSingleEleTrigger(triggerName, thresholdEle, caloId, caloIso, trkId, trkIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if (OpenHlt1ElectronPassed(thresholdEle[0], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
				       ) >= 1)	
	      
		triggerBit[it] = true;
	      }
	  }
      }

  else if (isSingleEleWPTrigger(triggerName, thresholdEle, caloId, caloIso, trkId, trkIso)){ 

    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
      { 
        if (prescaleResponse(menu, cfg, rcounter, it)) 
          { 
            if (OpenHlt1ElectronPassed(thresholdEle[0],  
                                       map_EGammaCaloId[caloId[0]], 
                                       map_EleCaloIso[caloIso[0]], 
                                       map_EleTrkId[trkId[0]], 
                                       map_EleTrkIso[trkIso[0]] 
                                       ) >= 1)   
               
	      triggerBit[it] = true; 
	  } 
      } 
  } 

  /*DoubleEle*/

  else if (isAsymDoubleEleTrigger(triggerName, thresholdEle, caloId, caloIso, trkId, trkIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if ((OpenHlt1ElectronPassed(thresholdEle[1], 
				       map_EGammaCaloId[caloId[1]],
				       map_EleCaloIso[caloIso[1]],
				       map_EleTrkId[trkId[1]],
				       map_EleTrkIso[trkIso[1]]
				       ) >= 2)	&&
		(OpenHlt1ElectronPassed(thresholdEle[0], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
				       ) >= 1)
		)
	      
		triggerBit[it] = true;
	      }
	  }
      }

  else if (isDoubleEleTrigger(triggerName, thresholdEle, caloId, caloIso, trkId, trkIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if (OpenHlt1ElectronPassed(thresholdEle[0], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
				       ) >= 2)	
	      
		triggerBit[it] = true;
	      }
	  }
      }


  /*SingleEle cross triggers*/

  else if (isSingleEleX_HTXTrigger(triggerName, thresholds, caloId, caloIso, trkId, trkIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if ((OpenHlt1ElectronPassed(thresholds[0], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
				       ) >= 1)	&&
		(OpenHltSumCorHTPassed(thresholds[1], 40.) == 1)
		)
	      
		triggerBit[it] = true;
	      }
	  }
      }



  /*SinglePhoton*/

  else if (isSinglePhotonTrigger(triggerName, thresholdPhoton, r9Id, caloId, photonIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if (OpenHlt1PhotonPassed(thresholdPhoton[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]
				     ) >= 1)	
	      
		triggerBit[it] = true;
	      }
	  }
      }

 /*DoublePhoton*/

//   else if (isDoublePhotonTrigger(triggerName, thresholdPhoton, r9Id, caloId, photonIso)){
    
//     if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
//       {
// 	if (prescaleResponse(menu, cfg, rcounter, it))
// 	  {
// 	    if (OpenHlt1PhotonPassed(thresholdPhoton[0],
// 				     map_PhotonR9ID[r9Id[0]],
// 				     map_EGammaCaloId[caloId[0]],
// 				     map_PhotonIso[photonIso[0]]
// 				     ) >= 2)	
	      
// 		triggerBit[it] = true;
// 	      }
// 	  }
//       }


//  version with EcalActiv , i.e. with no L1 seed

 else if (isDoublePhotonTrigger(triggerName, thresholdPhoton, r9Id, caloId, photonIso)){
   
   if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
     {
       std::vector<int> firstVector =  VecOpenHlt1PhotonPassed(thresholdPhoton[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]);
	      
	    if (firstVector.size()>=1){
	      std::vector<int> secondVector = VecOpenHlt1EcalActivPassed(thresholdPhoton[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]);

	      if (secondVector.size()>=1){
		for (unsigned int i=0; i<firstVector.size(); i++)
		  {
		    for (unsigned int j=0; j<secondVector.size() ; j++)
		      {
			if(abs(ohEcalActivEt[firstVector[i]] - ohPhotEt[secondVector[j]]) > 0.00001) 
			triggerBit[it] = true;
		      }
		  }
	      }
	    }
     }
 }

 /*DoublePhoton Asym*/

//   else if (isAsymDoublePhotonTrigger(triggerName, thresholdPhoton, r9Id, caloId, photonIso)){
    
//     if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
//       {
// 	if (prescaleResponse(menu, cfg, rcounter, it))
// 	  {
// 	     if ((OpenHlt1PhotonPassed(thresholdPhoton[1],
// 				     map_PhotonR9ID[r9Id[1]],
// 				     map_EGammaCaloId[caloId[1]],
// 				     map_PhotonIso[photonIso[1]]
// 				     ) >= 2) 	&&
	      
// 		(OpenHlt1PhotonPassed(thresholdPhoton[0],
// 				     map_PhotonR9ID[r9Id[0]],
// 				     map_EGammaCaloId[caloId[0]],
// 				     map_PhotonIso[photonIso[0]]
// 				      ) >= 1)
// 		 )

// 		triggerBit[it] = true;
// 	      }
// 	  }
//       }


//  version with EcalActiv , i.e. with no L1 seed

 else if (isAsymDoublePhotonTrigger(triggerName, thresholdPhoton, r9Id, caloId, photonIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {


	    std::vector<int> firstVector =  VecOpenHlt1PhotonPassed(thresholdPhoton[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]);
	      
	    if (firstVector.size()>=1){
	      std::vector<int> secondVector = VecOpenHlt1EcalActivPassed(thresholdPhoton[1],
				     map_PhotonR9ID[r9Id[1]],
				     map_EGammaCaloId[caloId[1]],
				     map_PhotonIso[photonIso[1]]);

	      if (secondVector.size()>=1){
		for (unsigned int i=0; i<firstVector.size(); i++)
		  {
		    for (unsigned int j=0; j<secondVector.size() ; j++)
		      {
			if(abs(ohEcalActivEt[firstVector[i]] - ohPhotEt[secondVector[j]] > 0.00001))
			triggerBit[it] = true;
		      }
		  }
	      }
	    }



	  }
      }
  }

  /* Single Jet */
  else if (triggerName.CompareTo("OpenHLT_L1SingleCenJet") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1SingleForJet") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1SingleTauJet") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1Jet6") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1Jet10") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      for (int i=0; i<NL1CenJet; i++)
		if (L1CenJetEt[i] >= 10.0)
		  rc++;
	      for (int i=0; i<NL1ForJet; i++)
		if (L1ForJetEt[i] >= 10.0)
		  rc++;
	      for (int i=0; i<NL1Tau; i++)
		if (L1TauEt [i] >= 10.0)
		  rc++;
	      if (rc > 0)
		triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L1Jet15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
	
  /*OpenHLT_JetX(U) paths*/
	
  else if (isJetXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1JetPassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isJetXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1CorJetPassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isJetX_NoJetIDTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1CorJetPassedNoJetID(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  
  /* DoubleJetX(U)_ForwardBackward */
	
  else if (isDoubleJetXU_ForwardBackwardTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int rc1 = 0;
	      int rc2 = 0;
				
	      // Loop over all oh jets, select events where both pT of a pair are above threshold and in HF+ and HF-
	      for (int i=0; i<NohJetCal; i++)
		{
		  if (ohJetCalPt[i]/0.7 > thresholds[0] && ohJetCalEta[i]
		      > 3.0 && ohJetCalEta[i] < 5.1)
		    { // Jet pT/eta cut
		      ++rc1;
		    }
		  if (ohJetCalPt[i]/0.7 > thresholds[0] && ohJetCalEta[i]
		      > -5.1 && ohJetCalEta[i] < -3.0)
		    { // Jet pT/eta cut
		      ++rc2;
		    }
		}
	      if (rc1!=0 && rc2!=0)
		rc=1;
	      if (rc > 0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDoubleJetX_ForwardBackwardTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int rc1 = 0;
	      int rc2 = 0;
				
	      // Loop over all oh jets, select events where both pT of a pair are above threshold and in HF+ and HF-
	      for (int i=0; i<NohJetCorCal; i++)
		{
		  if (ohJetCorCalPt[i] > thresholds[0]
		      && ohJetCorCalEta[i] > 3.0 && ohJetCorCalEta[i] < 5.1)
		    { // Jet pT/eta cut
		      ++rc1;
		    }
		  if (ohJetCorCalPt[i] > thresholds[0]
		      && ohJetCorCalEta[i] > -5.1 && ohJetCorCalEta[i]
		      < -3.0)
		    { // Jet pT/eta cut
		      ++rc2;
		    }
		}
	      if (rc1!=0 && rc2!=0)
		rc=1;
	      if (rc > 0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /* DiJetAveX(U) */
  else if (isDiJetAveXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltDiJetAvePassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDiJetAveXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltCorDiJetAvePassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
	
	
  /* Forward & MultiJet */
  else if (triggerName.CompareTo("OpenHLT_FwdJet20U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltFwdJetPassed(20.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /******ExlDiJetX(U)_HFAND**********/
  else if (isExclDiJetXU_HFANDTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      int rcDijetCand = 0;
	      double rcHFplusEnergy = 0;
	      double rcHFminusEnergy = 0;
				
	      // First loop over all jets and find a pair above threshold and with DeltaPhi/pi > 0.5
	      for (int i=0; i<NohJetCal; i++)
		{
		  if (ohJetCalPt[i]>thresholds[0])
		    { // Jet pT cut 
		      for (int j=0; j<NohJetCal && j!=i; j++)
			{
			  if (ohJetCalPt[j]>thresholds[0])
			    {
			      double Dphi=fabs(ohJetCalPhi[i]-ohJetCalPhi[j]);
			      if (Dphi>3.14159)
				Dphi=2.0*(3.14159)-Dphi;
			      if (Dphi>0.5*3.14159)
				{
				  rcDijetCand++;
				}
			    }
			}
		    }
		}
				
	      // Now ask for events with HF energy below threshold
	      if (rcDijetCand > 0)
		{
		  for (int i=0; i < NrecoTowCal; i++)
		    {
		      if ((recoTowEta[i] > 3.0) && (recoTowE[i] > 4.0))
			rcHFplusEnergy += recoTowE[i];
		      if ((recoTowEta[i] < -3.0) && (recoTowE[i] > 4.0))
			rcHFminusEnergy += recoTowE[i];
		    }
		}
				
	      // Now put them together
	      if ((rcDijetCand > 0) && (rcHFplusEnergy < 200) && (rcHFminusEnergy
								 < 200))
		triggerBit[it] = true;
	    }
	}
    }
	
  else if (isExclDiJetX_HFANDTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      int rcDijetCand = 0;
	      double rcHFplusEnergy = 0;
	      double rcHFminusEnergy = 0;
				
	      // First loop over all jets and find a pair above threshold and with DeltaPhi/pi > 0.5
	      for (int i=0; i<NohJetCorCal; i++)
		{
		  if (ohJetCorCalPt[i]>thresholds[0])
		    { // Jet pT cut 
		      for (int j=0; j<NohJetCorCal && j!=i; j++)
			{
			  if (ohJetCorCalPt[j]>thresholds[0])
			    {
			      double Dphi=fabs(ohJetCorCalPhi[i]-ohJetCorCalPhi[j]);
			      if (Dphi>3.14159)
				Dphi=2.0*(3.14159)-Dphi;
			      if (Dphi>0.5*3.14159)
				{
				  rcDijetCand++;
				}
			    }
			}
		    }
		}
				
	      // Now ask for events with HF energy below threshold
	      if (rcDijetCand > 0)
		{
		  for (int i=0; i < NrecoTowCal; i++)
		    {
		      if ((recoTowEta[i] > 3.0) && (recoTowE[i] > 4.0))
			rcHFplusEnergy += recoTowE[i];
		      if ((recoTowEta[i] < -3.0) && (recoTowE[i] > 4.0))
			rcHFminusEnergy += recoTowE[i];
		    }
		}
				
	      // Now put them together
	      if ((rcDijetCand > 0) && (rcHFplusEnergy < 200) && (rcHFminusEnergy
								 < 200))
		triggerBit[it] = true;
	    }
	}
    }
	
  /*******ExclDiJetX(U)_HFOR**********/
  else if (isExclDiJetXU_HFORTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      int rcDijetCand = 0;
	      double rcHFplusEnergy = 0;
	      double rcHFminusEnergy = 0;
				
	      // First loop over all jets and find a pair above threshold and with DeltaPhi/pi > 0.5 
	      for (int i=0; i<NohJetCal; i++)
		{
		  if (ohJetCalPt[i]>thresholds[0])
		    { // Jet pT cut  
		      for (int j=0; j<NohJetCal && j!=i; j++)
			{
			  if (ohJetCalPt[j]>thresholds[0])
			    {
			      double Dphi=fabs(ohJetCalPhi[i]-ohJetCalPhi[j]);
			      if (Dphi>3.14159)
				Dphi=2.0*(3.14159)-Dphi;
			      if (Dphi>0.5*3.14159)
				{
				  rcDijetCand++;
				}
			    }
			}
		    }
		}
				
	      // Now ask for events with HF energy below threshold 
	      if (rcDijetCand > 0)
		{
		  for (int i=0; i < NrecoTowCal; i++)
		    {
		      if ((recoTowEta[i] > 3.0) && (recoTowE[i] > 4.0))
			rcHFplusEnergy += recoTowE[i];
		      if ((recoTowEta[i] < -3.0) && (recoTowE[i] > 4.0))
			rcHFminusEnergy += recoTowE[i];
		    }
		}
				
	      // Now put them together 
	      if ((rcDijetCand > 0) && ((rcHFplusEnergy < 50) || (rcHFminusEnergy
								  < 50)))
		triggerBit[it] = true;
	    }
	}
    }
	
  else if (isExclDiJetX_HFORTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      int rcDijetCand = 0;
	      double rcHFplusEnergy = 0;
	      double rcHFminusEnergy = 0;
				
	      // First loop over all jets and find a pair above threshold and with DeltaPhi/pi > 0.5 
	      for (int i=0; i<NohJetCorCal; i++)
		{
		  if (ohJetCorCalPt[i]>thresholds[0])
		    { // Jet pT cut  
		      for (int j=0; j<NohJetCorCal && j!=i; j++)
			{
			  if (ohJetCorCalPt[j]>thresholds[0])
			    {
			      double Dphi=fabs(ohJetCorCalPhi[i]-ohJetCorCalPhi[j]);
			      if (Dphi>3.14159)
				Dphi=2.0*(3.14159)-Dphi;
			      if (Dphi>0.5*3.14159)
				{
				  rcDijetCand++;
				}
			    }
			}
		    }
		}
				
	      // Now ask for events with HF energy below threshold 
	      if (rcDijetCand > 0)
		{
		  for (int i=0; i < NrecoTowCal; i++)
		    {
		      if ((recoTowEta[i] > 3.0) && (recoTowE[i] > 4.0))
			rcHFplusEnergy += recoTowE[i];
		      if ((recoTowEta[i] < -3.0) && (recoTowE[i] > 4.0))
			rcHFminusEnergy += recoTowE[i];
		    }
		}
				
	      // Now put them together 
	      if ((rcDijetCand > 0) && ((rcHFplusEnergy < 50) || (rcHFminusEnergy
								  < 50)))
		triggerBit[it] = true;
	    }
	}
    }
  /********QuadJetX(U)********/
  else if (isQuadJetXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltQuadJetPassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isQuadJetXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltQuadCorJetPassed(thresholds[0])>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /**************************/
	
  else if (triggerName.CompareTo("OpenHLT_FwdJet40") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltFwdCorJetPassed(40.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_PentaJet25U20U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ((OpenHlt1JetPassed(25)>=4) && (OpenHlt1JetPassed(20)>=5))
		triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_QuadJet50_Jet40") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ((OpenHlt1CorJetPassed(50)>=4) && (OpenHlt1CorJetPassed(40)>=5))
		triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_QuadJet50_Jet40_Jet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ((OpenHlt1CorJetPassed(50)>=4) && (OpenHlt1CorJetPassed(40)>=5) && (OpenHlt1CorJetPassed(30)>=6))
		triggerBit[it] = true;
	    }
	}
    }
	
  /* MET, HT, SumHT, Razor, PT */
  else if (triggerName.CompareTo("OpenHLT_L1MET20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      triggerBit[it] = true;
	    }
	}
    }
	
  /***METX******/
  else if (isMETXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (recoMetCal > thresholds[0])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*******pfMHTX******/
	
  else if (ispfMHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (pfMHT > thresholds[0])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*******MeffX(U)********/
	
  else if (isMeffXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltMeffU(thresholds[0], 20.)==1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isMeffXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltMeff(thresholds[0], 40.)==1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /***********PTX(U)_X**************************/
	
  else if (isPT12U_XUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPT12U(thresholds[0], 50.)==1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isPT12_XTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPT12(thresholds[0], 50.)==1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*HTX and HTXU paths*/
	
  else if (isHTXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumHTPassed(thresholds[0], 20.) == 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(thresholds[0], 40.) == 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_SumET120") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (recoMetCalSum > 120.)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /****RX(U)_MRX(U)************/
  else if (isR0X_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[0], thresholds[1], 7, thresholds[2])>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
 else if (isR0X_MRX_BTagTrigger(triggerName, thresholds) )
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[0], thresholds[1], 7, thresholds[2])>0)
	    {
	      bool bjet=false;
	      for(int i=0; i<NohBJetL2Corrected;i++){
		if((ohBJetL2CorrectedEt[i] > thresholds[3] )
		   && (fabs(ohBJetL2CorrectedEta[i]) <  2.4)                                                                         
		   //		    && (ohBJetIPL25Tag[i] > 0.0)                                                                                    
		   && (ohBJetIPL3Tag[i]  > 6.0) ){
		  bjet = true;                                                                                                     
		}
	      }
	      if(bjet) 
		{
		  if (prescaleResponse(menu, cfg, rcounter, it)) 
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }

 else if (isRMRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRMRPassed(0.15,150,thresholds[0],-0.043,6., 7, 40.)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  else if (isDoublePhotonX_CaloIdL_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(0., thresholds[1], 7, 40.)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  if (OpenHlt1PhotonSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						    999.,
						    999., // Track iso barrel, Track iso endcap
						    999.,
						    999., // Track/pT iso barrel, Track/pT iso endcap
						    999.,
						    999., // H iso barrel, H iso endcap
						    999.,
						    999., // E iso barrel, E iso endcap
						    0.15,
						    0.10, // H/E barrel, H/E endcap
						    0.014,
						    0.035, // cluster shape barrel, cluster shape endcap
						    999.,//0.98,
						    999., // R9 barrel, R9 endcap
						    999.,
						    999., // Deta barrel, Deta endcap
						    999.,
						    999. // Dphi barrel, Dphi endcap
						    )>=2)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }

 else if (isDoublePhotonX_CaloIdL_R0X_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, thresholds[3])>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  if (OpenHlt1PhotonSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						    999.,
						    999., // Track iso barrel, Track iso endcap
						    999.,
						    999., // Track/pT iso barrel, Track/pT iso endcap
						    999.,
						    999., // H iso barrel, H iso endcap
						    999.,
						    999., // E iso barrel, E iso endcap
						    0.15,
						    0.10, // H/E barrel, H/E endcap
						    0.014,
						    0.035, // cluster shape barrel, cluster shape endcap
						    999.,//0.98,
						    999., // R9 barrel, R9 endcap
						    999.,
						    999., // Deta barrel, Deta endcap
						    999.,
						    999. // Dphi barrel, Dphi endcap
						    )>=2)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	

	
  else if (isPhotonX_R0X_MRXTrigger(triggerName, thresholds, r9Id, caloId, photonIso))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, thresholds[3])>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  if (OpenHlt1PhotonPassed(thresholds[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]
				     ) >= 1)		
		    {
		      triggerBit[it] = true;
		       }
		}
	    }
	}
    }
	
  else if (isDoublePhotonX_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(0., thresholds[1], 7, 40.)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  if (OpenHlt1PhotonSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						    999.,
						    999., // Track iso barrel, Track iso endcap
						    999.,
						    999., // Track/pT iso barrel, Track/pT iso endcap
						    999.,
						    999., // H iso barrel, H iso endcap
						    999.,
						    999., // E iso barrel, E iso endcap
						    0.15,
						    0.10, // H/E barrel, H/E endcap
						    999.,
						    999., // cluster shape barrel, cluster shape endcap
						    999.,//0.98,
						    999., // R9 barrel, R9 endcap
						    999.,
						    999., // Deta barrel, Deta endcap
						    999.,
						    999. // Dphi barrel, Dphi endcap
						    )>=2)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  else if (isDoublePhotonX_R0X_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, thresholds[3])>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  if (OpenHlt1PhotonSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						    999.,
						    999., // Track iso barrel, Track iso endcap
						    999.,
						    999., // Track/pT iso barrel, Track/pT iso endcap
						    999.,
						    999., // H iso barrel, H iso endcap
						    999.,
						    999., // E iso barrel, E iso endcap
						    0.15,
						    0.10, // H/E barrel, H/E endcap
						    999.,
						    999., // cluster shape barrel, cluster shape endcap
						    999.,//0.98,
						    999., // R9 barrel, R9 endcap
						    999.,
						    999., // Deta barrel, Deta endcap
						    999.,
						    999. // Dphi barrel, Dphi endcap
						    )>=2)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
  else if (isMuX_R0X_MRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 0)>=1)
		{
	     	  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, 40.)>0)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (isMuX_MR0XTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 0)>=1)
		{
		  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, 40.)>0)
		    triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isEleX_R0X_MRXTrigger(triggerName, thresholds, caloId, caloIso, trkId, trkIso))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ( OpenHlt1ElectronPassed(thresholds[0], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
				       ) >= 1)	
		{
		  if (OpenHltRPassed(thresholds[1], thresholds[2], 7, 40.)>0)
		    triggerBit[it] = true;
		}
	    }
	}
    }


  else if (isEleLX_MR0XTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ( (OpenHlt1ElectronSamHarperPassed(thresholds[1], 0, // ET, L1isolation 
						    999.,
						    999., // Track iso barrel, Track iso endcap 
						    0.2,
						    0.2, // Track/pT iso barrel, Track/pT iso endcap 
						    0.2,
						    0.2, // H/ET iso barrel, H/ET iso endcap 
						    0.2,
						    0.2, // E/ET iso barrel, E/ET iso endcap 
						    0.15,
						    0.10, // H/E barrel, H/E endcap 
						    0.014,
						    0.035, // cluster shape barrel, cluster shape endcap 
						    999.,//0.98,
						    999.,//1.0, // R9 barrel, R9 endcap 
						    0.01,
						    0.01, // Deta barrel, Deta endcap 
						    0.15,
						    0.10 // Dphi barrel, Dphi endcap 
						    )>=1))
		{
		  if (OpenHltRPassed(0., thresholds[0], 7, 40.)>0)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (isEleTX_MR0XTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(thresholds[2], 0, // ET, L1isolation                                                                                                                                                                           
						  999.,
						  999., // Track iso barrel, Track iso endcap     
						  0.2,
						  0.2, // Track/pT iso barrel, Track/pT iso endcap
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap
						  0.1,
						  0.075, // H/E barrel, H/E endcap                
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap                        
						  0.008,
						  0.008, // Deta barrel, Deta endcap                  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap
						  )>=1)
		{
		  if (OpenHltRPassed(0., thresholds[0], 7, 40.)>0)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (isR0XTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(thresholds[0], 0., 7, 40.)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isMRXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHltRPassed(0., thresholds[0], 7, 40.)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	

  /******FatJets********/
//   FATJETS
   else if(isFatJetMass(triggerName, thresholds))
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if(OpenHltFatJetPassed(30., thresholds[1], thresholds[2], thresholds[0])) 
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

/// FATJETS BTAG
   else if(isFatJetMassBTag(triggerName, thresholds))
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if(OpenHltFatJetPassed(30., thresholds[1], thresholds[2], thresholds[0])) 
		 {
		   int rc = 0;
		   int max = (NohBJetL2Corrected > 2) ? 2 : NohBJetL2Corrected;
		   for (int i = 0; i < max; i++)
		     {
		       if (ohBJetL2CorrectedEt[i] > thresholds[0])
			 { // ET cut 
			   if (ohBJetIPL25Tag[i] > 2.5)
			     { // Level 2.5 b tag 
			       if (ohBJetIPL3Tag[i] > 3.5)
				 { // Level 3 b tag 
				   rc++;
				 }
			     }
			 }
		     }
		   if (rc >= 1)
		     {
		       triggerBit[it] = true;
		     }
		 }
	     }
	 }
     }
/// SingleJet BTAG
   else if (isJetX_BTagIPTrigger(triggerName, thresholds))
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1CorJetPassed(thresholds[0])>=1)
            {
	      int rc = 0;
	      int max = (NohBJetL2Corrected > 2) ? 2 : NohBJetL2Corrected;
	      for (int i = 0; i < max; i++)
		{
                  if (ohBJetL2CorrectedEt[i] > thresholds[0])
		    { // ET cut 
		      if (ohBJetIPL25Tag[i] > 2.5)
			{ // Level 2.5 b tag 
			  if (ohBJetIPL3Tag[i] > 3.5)
			    { // Level 3 b tag 
			      rc++;
			    }
			}
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	 }
      }
   }




	
  /****BTagIP_HTX************/
	
  else if (isHTX_CentralJetX_BTagIPTrigger(triggerName, thresholds)) {
    //if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second == 1){
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0){
			
      if (prescaleResponse(menu, cfg, rcounter, it) ) {
	if (OpenHltSumCorHTPassed(thresholds[0], 40.) == 1){
					
	  int rc = 0;
	  int max = (NohBJetL2Corrected > 6) ? 6 : NohBJetL2Corrected;
	  for (int i = 0; i < max; i++){
	    if (ohBJetL2CorrectedEt[i] > thresholds[1] && fabs(ohBJetL2CorrectedEta[i]) < 3.0)
	      { // ET cut 
		//if (ohBJetIPL25Tag[i] > 0.0)
		//{ // Level 2.5 b tag 
		if (ohBJetIPL3Tag[i] > 4.0)
		  { // Level 3 b tag 
		    rc++;
		  }
		//}
	      }
	  }
	  if (rc >= 1){
	    triggerBit[it] = true;
						
	  }	   
	}	       
      }
    }
  }
	
  /****BTagIP_pfMHTX_HTX************/
	
  else if (isBTagIP_pfMHTX_HTXTrigger(menu->GetTriggerName(it), thresholds)){
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu, cfg, rcounter, it)){
	if (OpenHltSumCorHTPassed(thresholds[0], 40)>=1 && pfMHT >=thresholds[2]){
					
	  int rc = 0;
	  int max = (NohBJetL2Corrected > 6) ? 6 : NohBJetL2Corrected;
	  for (int i = 0; i < max; i++){
	    if (ohBJetL2CorrectedEt[i] > thresholds[1] && fabs(ohBJetL2CorrectedEta[i]) < 3.0)
	      { // ET cut 
		//if (ohBJetIPL25Tag[i] > 0.0)
		//{ // Level 2.5 b tag 
		if (ohBJetIPL3Tag[i] > 4.0)
		  { // Level 3 b tag 
		    rc++;
		  }
		// }
	      }
	  }
	  if (rc >= 1){
	    triggerBit[it] = true;
	  }	   
	}
      }
    }
  }
	
  /****pfMHTX_HTX************/
	
  else if (ispfMHTX_HTXTrigger(menu->GetTriggerName(it), thresholds)){
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu, cfg, rcounter, it)){
	if (OpenHltSumCorHTPassed(thresholds[0], 40)>=1 && pfMHT >=thresholds[1]){
	  triggerBit[it] = true;
	}
      }
    }
  }
	
  /* Muons */
  else if (isL1SingleMuXTrigger(triggerName))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (true)
		{ // passthrough
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isL1DoubleMuXTrigger(triggerName))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (true)
		{ // passthrough
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  else if (isL2SingleMuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonPassed(0.0, thresholds[0], 9999.0) > 0)
		triggerBit[it] = true;
	    }
	}
    }
	
  else if (isL2DoubleMuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonPassed(0.0, thresholds[0], 9999.0) > 1)
		triggerBit[it] = true;
	    }
	}
    }
  else if (isMuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }	
  else if (isMuX_eta2pXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 0, thresholds[1], thresholds[1])>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (isIsoMuX_eta2pXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 1, thresholds[1], thresholds[1])>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
 else if (isIsoMuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 1)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }	
 
  else if (triggerName.CompareTo("OpenHLT_Mu0_L1MuOpen") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., 0., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu0_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., 0., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu3_L1MuOpen") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_L1MuOpen") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
 
  else if (triggerName.CompareTo("OpenHLT_Mu3Mu0") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., 0., 2., 0)>=2 && OpenHlt1MuonPassed(
										 0.,
										 3.,
										 3.,
										 2.,
										 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	


  else if (triggerName.CompareTo("OpenHLT_L2Mu0_NoVertex") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonNoVertexPassed(0., 0., 9999.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L2DoubleMu45_NoVertex") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
              if (OpenHlt1L2MuonNoVertexPassed(0., 45., 9999.)>=2) 
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_L2DoubleMu23_NoVertex") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonNoVertexPassed(0., 23., 9999.)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L2DoubleMu35_NoVertex_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonNoVertexPassed(0., 35., 9999.)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L2DoubleMu30_NoVertex") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1L2MuonNoVertexPassed(0., 30., 9999.)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu50_NoVertex") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 15., 50., 9999., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (isDoubleMuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., thresholds[0], 2., 0)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  else if (triggerName.CompareTo("OpenHLT_DoubleMu6_Acoplanarity03") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      float ptl2 = 3.5;
	      float ptl3 = 6.0;
	      float drl3 = 2.0;
	      float etal3 = 2.15;
	      float etal2 = 2.15;
	      float deltaphil3 = 0.3;
	      float deltaptl3 = 3.0;
				
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
		      if (fabs(i[ohMuL3Eta]) < etal3 && fabs(ohMuL3Eta[j]) < etal3)
			{ // L3 eta cut   
			  if (ohMuL3Pt[i] > ptl3 && ohMuL3Pt[j] > ptl3)
			    { // L3 pT cut         
			      if (ohMuL3Dr[i] < drl3 && ohMuL3Dr[j] < drl3)
				{ // L3 DR cut 
				  if ((ohMuL3Chg[i] * ohMuL3Chg[j]) < 0)
				    { // opposite charge
				      float deltaphi = fabs(ohMuL3Phi[i]-ohMuL3Phi[j]);
				      if (deltaphi > 3.14159)
					deltaphi = (2.0 * 3.14159) - deltaphi;
										
				      deltaphi = 3.14159 - deltaphi;
				      if (deltaphi < deltaphil3)
					{
					  float deltapt = fabs(ohMuL3Pt[i]-ohMuL3Pt[j]);
					  if(deltapt < deltaptl3)
					    {
					      int l2match1 = ohMuL3L2idx[i];
					      int l2match2 = ohMuL3L2idx[j];
					      
					      if ( (fabs(ohMuL2Eta[l2match1]) < etal2)
						   && (fabs(ohMuL2Eta[l2match2]) < etal2))
						{ // L2 eta cut  
						  if ( (ohMuL2Pt[l2match1] > ptl2)
						       && (ohMuL2Pt[l2match2] > ptl2))
						    { // L2 pT cut  
						      rc++;
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	      if (rc >=1)
		triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu7_Acoplanarity03") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
              int rc = 0; 
              float ptl2 = 3.5; 
              float ptl3 = 7.0; 
              float drl3 = 2.0; 
              float etal3 = 2.15; 
              float etal2 = 2.15; 
              float deltaphil3 = 0.3; 
              float deltaptl3 = 3.0; 
              for (int i=0; i<NohMuL3; i++) 
                { 
                  for (int j=i+1; j<NohMuL3; j++) 
                    { 
                      if (fabs(i[ohMuL3Eta]) < etal3 && fabs(ohMuL3Eta[j]) < etal3) 
                        { // L3 eta cut    
                          if (ohMuL3Pt[i] > ptl3 && ohMuL3Pt[j] > ptl3) 
                            { // L3 pT cut          
                              if (ohMuL3Dr[i] < drl3 && ohMuL3Dr[j] < drl3) 
                                { // L3 DR cut  
                                  if ((ohMuL3Chg[i] * ohMuL3Chg[j]) < 0) 
                                    { // opposite charge 
                                      float deltaphi = fabs(ohMuL3Phi[i]-ohMuL3Phi[j]); 
                                      if (deltaphi > 3.14159) 
                                        deltaphi = (2.0 * 3.14159) - deltaphi; 
                                                                                 
                                      deltaphi = 3.14159 - deltaphi; 
                                      if (deltaphi < deltaphil3) 
                                        { 
                                          float deltapt = fabs(ohMuL3Pt[i]-ohMuL3Pt[j]); 
                                          if(deltapt < deltaptl3) 
                                            { 
                                              int l2match1 = ohMuL3L2idx[i]; 
                                              int l2match2 = ohMuL3L2idx[j]; 
                                               
                                              if ( (fabs(ohMuL2Eta[l2match1]) < etal2) 
                                                   && (fabs(ohMuL2Eta[l2match2]) < etal2)) 
                                                { // L2 eta cut   
                                                  if ( (ohMuL2Pt[l2match1] > ptl2) 
                                                       && (ohMuL2Pt[l2match2] > ptl2)) 
                                                    { // L2 pT cut   
                                                      rc++; 
                                                    } 
                                                } 
                                            } 
                                        } 
                                    } 
                                } 
                            } 
                        } 
                    } 
                } 
              if (rc >=1) 
                triggerBit[it] = true; 
            } 
        } 
    } 

  else if (isMuX_MuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
	      if (OpenHlt1MuonPassed(0., 0., thresholds[3], 2., 0)>=2 && OpenHlt1MuonPassed(thresholds[0], thresholds[1], thresholds[2], 2., 0)>=1)
		{
                  triggerBit[it] = true;
                }
            }
        }
    }
	
  else if (triggerName.CompareTo("OpenHLT_L1DoubleMuOpen") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (1)
		{ // Pass through
		  triggerBit[it] = true;
		}
	    }
	}
    }
 
  //non standard muon paths
  else if (triggerName.CompareTo("OpenHLT_TripleMu5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu5_IsoMu5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=3 && OpenHlt1MuonPassed(
										 3.,
										 3.,
										 5.,
										 2.,
										 1)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_TripleMu5_2IsoMu5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=3 && OpenHlt1MuonPassed(
										 3.,
										 4.,
										 5.,
										 2.,
										 1)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_TripleIsoMu5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 1)>=3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_TripleMu_2IsoMu5_1Mu8") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 1)>=3&&OpenHlt1MuonPassed(
									       3.,
									       4.,
									       8.,
									       2.,
									       0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
 
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
              if (OpenHlt1MuonPassed(0., 10., 17., 2., 1, 2.1, 2.1)>=1) 
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
 
  else if (triggerName.CompareTo("OpenHLT_IsoMu20_eta2p1") == 0)  
    {  
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)  
        {  
          if (prescaleResponse(menu, cfg, rcounter, it))  
            {  
              if (OpenHlt1MuonPassed(0., 12., 20., 2., 1, 2.1, 2.1)>=1)  
                {  
                  triggerBit[it] = true;  
                }  
            }  
        }  
    }  
	
  /* Quarkonia */
  else if (triggerName.CompareTo("OpenHLT_DoubleMu2_Bs") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassPassed(0., 0., 2., 2., 0, 4.8, 6.0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu3_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassPassed(0., 0., 3., 2., 0, 2.5, 4.0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu3_Quarkonium") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassPassed(0., 0., 3., 2., 0, 1.5, 14.0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Onia") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[2] = { 2.6, 7.5 };
	      double massMaxPix[2] = { 3.6, 12.0 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[2] = { 2.8, 8.5 };
	      double massMaxTrack[2] = { 3.4, 11.0 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed(
				       ptPix,
				       pPix,
				       etaPix,
				       DxyPix,
				       DzPix,
				       NHitsPix,
				       normChi2Pix,
				       massMinPix,
				       massMaxPix,
				       DzMuonPix,
				       checkChargePix) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed(
				       ptTrack,
				       pTrack,
				       etaTrack,
				       DxyTrack,
				       DzTrack,
				       NHitsTrack,
				       normChi2Track,
				       massMinTrack,
				       massMaxTrack,
				       DzMuonTrack,
				       checkChargeTrack))
		{ //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	      //        if (GetIntRandom() % menu->GetPrescale(it) == 0) { triggerBit[it] = true; }
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu0_Track0_Ups") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 7.5 };
	      double massMaxPix[1] = { 12.0 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 8.5 };
	      double massMaxTrack[1] = { 11.0 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 0., 0., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_Ups(
					   ptPix,
					   pPix,
					   etaPix,
					   DxyPix,
					   DzPix,
					   NHitsPix,
					   normChi2Pix,
					   massMinPix,
					   massMaxPix,
					   DzMuonPix,
					   checkChargePix,
					   7) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_Ups(
					   ptTrack,
					   pTrack,
					   etaTrack,
					   DxyTrack,
					   DzTrack,
					   NHitsTrack,
					   normChi2Track,
					   massMinTrack,
					   massMaxTrack,
					   DzMuonTrack,
					   checkChargeTrack,
					   7))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu3_Track0_Ups") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 7.5 };
	      double massMaxPix[1] = { 12.0 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 8.5 };
	      double massMaxTrack[1] = { 11.0 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_Ups(
					   ptPix,
					   pPix,
					   etaPix,
					   DxyPix,
					   DzPix,
					   NHitsPix,
					   normChi2Pix,
					   massMinPix,
					   massMaxPix,
					   DzMuonPix,
					   checkChargePix,
					   8) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_Ups(
					   ptTrack,
					   pTrack,
					   etaTrack,
					   DxyTrack,
					   DzTrack,
					   NHitsTrack,
					   normChi2Track,
					   massMinTrack,
					   massMaxTrack,
					   DzMuonTrack,
					   checkChargeTrack,
					   8))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_Track0_Ups") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 7.5 };
	      double massMaxPix[1] = { 12.0 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 8.5 };
	      double massMaxTrack[1] = { 11.0 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_Ups(
					   ptPix,
					   pPix,
					   etaPix,
					   DxyPix,
					   DzPix,
					   NHitsPix,
					   normChi2Pix,
					   massMinPix,
					   massMaxPix,
					   DzMuonPix,
					   checkChargePix,
					   9) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_Ups(
					   ptTrack,
					   pTrack,
					   etaTrack,
					   DxyTrack,
					   DzTrack,
					   NHitsTrack,
					   normChi2Track,
					   massMinTrack,
					   massMaxTrack,
					   DzMuonTrack,
					   checkChargeTrack,
					   9))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu0_Track0_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 0., 0., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    0) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    0))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu3_Track0_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    5) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    5))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu3_Track3_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 3.;
	      double pTrack = 5.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    5) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    5))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  else if (triggerName.CompareTo("OpenHLT_Mu3_Track5_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 5.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(0., 3., 3., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    5) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    5))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_Track0_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 0.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    6) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    6))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_Track5_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 5.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    6) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    6))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu7_Track5_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 5.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(3., 4., 7., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    6) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    6))
		{
		  //check the L3Mu + tracker track
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu7_Track7_Jpsi") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      //       cout << "checking for Onia " << endl;
	      //variables for pixel cuts
	      double ptPix = 0.;
	      double pPix = 3.;
	      double etaPix = 999.;
	      double DxyPix = 999.;
	      double DzPix = 999.;
	      int NHitsPix = 3;
	      double normChi2Pix = 999999999.;
	      double massMinPix[1] = { 2.6 };
	      double massMaxPix[1] = { 3.6 };
	      double DzMuonPix = 999.;
	      bool checkChargePix = false;
	      //variables for tracker track cuts
	      double ptTrack = 7.;
	      double pTrack = 3.;
	      double etaTrack = 999.;
	      double DxyTrack = 999.;
	      double DzTrack = 999.;
	      int NHitsTrack = 5;
	      double normChi2Track = 999999999.;
	      double massMinTrack[1] = { 2.8 };
	      double massMaxTrack[1] = { 3.4 };
	      double DzMuonTrack = 0.5;
	      bool checkChargeTrack = true;
	      if ((OpenHlt1MuonPassed(3., 4., 7., 2., 0)>=1) && //check the L3 muon
		  OpenHltMuPixelPassed_JPsi(
					    ptPix,
					    pPix,
					    etaPix,
					    DxyPix,
					    DzPix,
					    NHitsPix,
					    normChi2Pix,
					    massMinPix,
					    massMaxPix,
					    DzMuonPix,
					    checkChargePix,
					    6) && //check the L3Mu + pixel
		  OpenHltMuTrackPassed_JPsi(
					    ptTrack,
					    pTrack,
					    etaTrack,
					    DxyTrack,
					    DzTrack,
					    NHitsTrack,
					    normChi2Track,
					    massMinTrack,
					    massMaxTrack,
					    DzMuonTrack,
					    checkChargeTrack,
					    6))
		{
		  {
		    //check the L3Mu + tracker track
		    triggerBit[it] = true;
		  }
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_DoubleMu0_Quarkonium") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=0; j<NohMuL3 && j != i; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      if (diMuMass > 2.5 && diMuMass < 14.5 && dimuCharge == 0)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_DoubleMu0_Quarkonium_LS") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=0; j<NohMuL3 && j != i; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      if (diMuMass > 2.5 && diMuMass < 14.5 && dimuCharge != 0)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  /*dimuon for Higgs*/
  else if (triggerName.CompareTo("OpenHLT_DiMuon6_LowMass") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxHiggsPassed(3., 3., 4., 2., 1, 5., 14.5, 2.5, 6., 0.5, 10, -2, 6.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
       }
    }
  		
	
  /*MuOnia 1e33*/
  else if (triggerName.CompareTo("OpenHLT_Dimuon6p5_Barrel_PsiPrime_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      float diMuPt = diMu.Pt();
		      float diMuRapidity = diMu.Rapidity();
		      if (diMuMass > 3.35 && diMuMass < 4.0 && dimuCharge == 0 && diMuPt > 6.5 && fabs(diMuRapidity) < 1.3)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Dimuon6p5_Jpsi_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      float diMuPt = diMu.Pt();
		      if (diMuMass > 2.5 && diMuMass < 4.0 && dimuCharge == 0 && diMuPt > 6.5)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }   
	
  else if (triggerName.CompareTo("OpenHLT_Dimuon6p5_Barrel_Jpsi_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      float diMuPt = diMu.Pt();
		      float diMuRapidity = diMu.Rapidity();
		      if (diMuMass > 2.5 && diMuMass < 4.0 && dimuCharge == 0 && diMuPt > 6.5 && fabs(diMuRapidity) < 1.3)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }    
	
  else if (triggerName.CompareTo("OpenHLT_Dimuon0_Barrel_Upsilon_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      float diMuRapidity = diMu.Rapidity();
		      if (diMuMass > 8.5 && diMuMass < 11.5 && dimuCharge == 0 && fabs(diMuRapidity) < 1.3)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }    
	
	
	
  else if (triggerName.CompareTo("OpenHLT_Dimuon6p5_LowMass_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      TLorentzVector mu1;
	      TLorentzVector mu2;
	      TLorentzVector diMu;
	      const double muMass = 0.105658367;
	      int rc = 0;
	      for (int i=0; i<NohMuL3; i++)
		{
		  for (int j=i+1; j<NohMuL3; j++)
		    {
						
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      int dimuCharge = (int) (ohMuL3Chg[i] + ohMuL3Chg[j]);
		      float diMuMass = diMu.M();
		      float diMuPt = diMu.Pt();
		      if (diMuMass > 1.0 && diMuMass < 5.0 && dimuCharge == 0 && diMuPt > 6.5)
			rc++;
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }    

  /*MuOnia 1e33 and 2e33*/
 else if (triggerName.CompareTo("OpenHLT_Dimuon0_Jpsi_Muon_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt3MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.8, 3.35, 2.5, 0.0, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon0_Jpsi_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.8, 3.35, 2.5, 0.0, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon0_Upsilon_Muon_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt3MuonOSMassVtxPassed(0., 0., 0., 2., 0, 8.5, 11.5, 2.5, 0.0, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon0_Upsilon_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 8.5, 11.5, 2.5, 0.0, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon10_Jpsi_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.8, 3.35, 2.5, 9.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon15_Jpsi_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.8, 3.35, 2.5, 14.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon4_Bs_Barrel_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassDCAPassed(0., 0., 2., 2., 0, 4.8, 6.0, 1.5, 3.9, 0.5)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon5_Upsilon_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 8.5, 11.5, 2.5, 4.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon10_Upsilon_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 8.5, 11.5, 2.5, 9.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon6_Bs_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassDCAPassed(0., 0., 2., 2., 0, 4.8, 6.0, 2.5, 5.9, 0.5)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon7_Jpsi_Displaced_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.9, 3.3, 2.5, 6.9, 0.5, 7.8794, 0.9, 3.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu3p5_Jpsi_Displaced_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 3.5, 2., 0, 2.9, 3.3, 2.2, 6.9, 0.5, 2.7055, 0.9, 3.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon7_Jpsi_X_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.95, 3.25, 2.5, 6.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon12_Jpsi_X_Barrel_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 2.95, 3.25, 2.5, 11.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon7_LowMass_Displaced_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 3., 2., 0, 1.0, 4.8, 2.2, 6.9, 0.5, 3.8414, 0.9, 3.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_DoubleMu4_LowMass_Displaced_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 4., 2., 0, 1.0, 4.8, 2.2, 6.9, 0.5, 2.0722, 0.9, 3.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon7_PsiPrime_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 3.35, 4.05, 2.5, 6.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Dimuon12_PsiPrime_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2MuonOSMassVtxPassed(0., 0., 0., 2., 0, 3.35, 4.05, 2.5, 11.9, 0.5, 7.8794, -2, 0.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
		
	
  /* Electrons */
	
  else if (isL1SingleEGXTrigger(triggerName))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (true)
		{ // passthrough      
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_L1DoubleEG5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (true)
		{ // passthrough        
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele8_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  999.,
						  999., // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  999,
						  999, // Deta barrel, Deta endcap  
						  999,
						  999 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdL_CaloIsoVL_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999,
						  999, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap  
						  0.98,
						  1.0, // R9 barrel, R9 endcap  
						  999,
						  999, // Deta barrel, Deta endcap  
						  999,
						  999 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-05-11: promoted to v4. R9 removed from prev. version. 
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdL_TrkIdVL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  0.01,
						  0.01, // Deta barrel, Deta endcap  
						  0.15,
						  0.1 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

	//2011-05-22:
else if (triggerName.CompareTo("OpenHLT_DoubleEle8_CaloIdT_TrkIdT_v1") == 0)//new
{
	if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
		if (prescaleResponse(menu, cfg, rcounter, it))
		{
            if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
												999.,
												999., // Track iso barrel, Track iso endcap  
												999,
												999, // Track/pT iso barrel, Track/pT iso endcap  
												999.,
												999., // H/ET iso barrel, H/ET iso endcap  
												999.,
												999., // E/ET iso barrel, E/ET iso endcap  
												0.1,
												0.075, // H/E barrel, H/E endcap 
												0.011,
												0.031, // cluster shape barrel, cluster shape endcap 
												999.,
												999., // R9 barrel, R9 endcap 
												0.008,
												0.008, // Deta barrel, Deta endcap
												0.07,
												0.05 // Dphi barrel, Dphi endcap 
												)>=2 )
            {
				triggerBit[it] = true;
            }
		}
	}
}

	
  //JH
 else if (triggerName.CompareTo("OpenHLT_Ele22_CaloIdL_CaloIsoVL_Ele15_HFT") == 0) {
   if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1) {
     if (prescaleResponse(menu, cfg, rcounter, it)) {
       // 17 GeV Electron
       if (
	   OpenHlt1ElectronSamHarperPassed(
					   22., 0, // ET, L1isolation
					   999., 999., // Track iso barrel, Track iso endcap
					   999., 999., // Track/pT iso barrel, Track/pT iso endcap
					   .2, .2, // H/ET iso barrel, H/ET iso endcap
					   .2, .2, // E/ET iso barrel, E/ET iso endcap
					   0.15, 0.1,  // H/E barrel, H/E endcap
					   0.014, 0.035, // cluster shape barrel, cluster shape endcap
					   999, 999,  // R9 barrel, R9 endcap
					   999, 999,   // Deta barrel, Deta endcap
					   999, 999    // Dphi barrel, Dphi endcap
					   ) >=1
	   ) {
	 // HF e code
	 float MinPt  = 15.0;
	 float MaxEta = 5.0;
	 float MinEta = -5.0;
	 float e9e25Min = 0.92; // e9e25 > 0.92
	 float e1e9Min = 0.6; // e1e9 > 0.6
	 float Min2DCut = 0.2; // 2DCut > 0.2

	 for(int i = 0; i < NohHFECALClus; i++){ // These arrays are kept in sync when creating openhlt.root
	   if (
                               ohHFEleClustere9e25[i] > e9e25Min &&
                               ohHFEleClustere1e9[i]  > e1e9Min  &&
                               ohHFEleCluster2Dcut[i] > Min2DCut &&
                               ohHFEleEta[i] < MaxEta &&
                               ohHFEleEta[i] > MinEta &&
                               ohHFElePt[i]  > MinPt
			       ) {
	     triggerBit[it] = true; // Fired
	   }
	 }
       }
     }
   }
 }

  //end JH	

  // 2011-05-11: promoted to v4. Removed R9 cut.
  else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdL_CaloIsoVL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(17., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999,
						  999, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap 
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap 
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  999,
						  999, // Deta barrel, Deta endcap  
						  999,
						  999 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(15., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.125,
						  0.075, // Track/pT iso barrel, Track/pT iso endcap  
						  0.125,
						  0.075, // H/ET iso barrel, H/ET iso endcap  
						  0.125,
						  0.075, // E/ET iso barrel, E/ET iso endcap  
						  0.05,
						  0.05, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele52_CaloIdVT_TrkIdT_v1")  //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(52., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.05,
						  0.05, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele27_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1") //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(27., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.2,
						  0.2, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  0.024,
						  0.040, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  0.01,
						  0.01, // Deta barrel, Deta endcap  
						  0.15,
						  0.1 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  //2011-06-08
  else if (triggerName.CompareTo("OpenHLT_Ele100_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1") 
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1ElectronSamHarperPassed(100., 0, // ET, L1isolation  
                                                  999.,
                                                  999., // Track iso barrel, Track iso endcap  
                                                  0.2,
                                                  0.2, // Track/pT iso barrel, Track/pT iso endcap  
                                                  0.2,
                                                  0.2, // H/ET iso barrel, H/ET iso endcap  
                                                  0.2,
                                                  0.2, // E/ET iso barrel, E/ET iso endcap  
                                                  0.15,
                                                  0.10, // H/E barrel, H/E endcap  
                                                  0.024,
                                                  0.040, // cluster shape barrel, cluster shape endcap  
                                                  999.,
                                                  999., // R9 barrel, R9 endcap  
                                                  0.01,
                                                  0.01, // Deta barrel, Deta endcap  
                                                  0.15,
                                                  0.1 // Dphi barrel, Dphi endcap  
                                                  )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }

	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1")  //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(32., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.2,
						  0.2, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  0.024,
						  0.040, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  0.010,
						  0.010, // Deta barrel, Deta endcap  
						  0.15,
						  0.10 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // JH
  else if (triggerName.CompareTo("OpenHLT_Ele42_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL")
	   == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
              if (OpenHlt1ElectronSamHarperPassed(42., 0, // ET, L1isolation   
                                                  999., 
                                                  999., // Track iso barrel, Track iso endcap   
                                                  0.2, 
                                                  0.2, // Track/pT iso barrel, Track/pT iso endcap   
                                                  0.2, 
                                                  0.2, // H/ET iso barrel, H/ET iso endcap   
                                                  0.2, 
                                                  0.2, // E/ET iso barrel, E/ET iso endcap   
                                                  0.15, 
                                                  0.10, // H/E barrel, H/E endcap   
                                                  0.024, 
                                                  0.040, // cluster shape barrel, cluster shape endcap   
                                                  999., 
                                                  999., // R9 barrel, R9 endcap   
                                                  0.01, 
                                                  0.01, // Deta barrel, Deta endcap   
                                                  0.15, 
                                                  0.1 // Dphi barrel, Dphi endcap   
                                                  )>=1) 
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT")
           == 0)  
    {  
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)  
        {  
          if (prescaleResponse(menu, cfg, rcounter, it))  
            {  
              if (OpenHlt1ElectronSamHarperPassed(42., 0, // ET, L1isolation    
                                                  999., 
                                                  999., // Track iso barrel, Track iso endcap   
                                                  0.125, 
                                                  0.075, // Track/pT iso barrel, Track/pT iso endcap   
                                                  0.125, 
                                                  0.075, // H/ET iso barrel, H/ET iso endcap   
                                                  0.125, 
                                                  0.075, // E/ET iso barrel, E/ET iso endcap   
                                                  0.05, 
                                                  0.05, // H/E barrel, H/E endcap   
                                                  0.011, 
                                                  0.031, // cluster shape barrel, cluster shape endcap   
                                                  999.,//0.98, 
                                                  999.,//1.0, // R9 barrel, R9 endcap   
                                                  0.008, 
                                                  0.008, // Deta barrel, Deta endcap   
                                                  0.07, 
                                                  0.05 // Dphi barrel, Dphi endcap   
                                                  )>=1)  
                {  
                  triggerBit[it] = true;  
                }  
            }  
        }  
    }  
  else if (triggerName.CompareTo("OpenHLT_Ele52_CaloIdVT_TrkIdT")
           == 0)   
    {   
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)   
        {   
          if (prescaleResponse(menu, cfg, rcounter, it))   
            {   
              if (OpenHlt1ElectronSamHarperPassed(52., 0, // ET, L1isolation     
                                                  999., 
                                                  999., // Track iso barrel, Track iso endcap   
                                                  999., 
                                                  999., // Track/pT iso barrel, Track/pT iso endcap   
                                                  999., 
                                                  999., // H/ET iso barrel, H/ET iso endcap   
                                                  999., 
                                                  999., // E/ET iso barrel, E/ET iso endcap   
                                                  0.05, 
                                                  0.05, // H/E barrel, H/E endcap   
                                                  0.011, 
                                                  0.031, // cluster shape barrel, cluster shape endcap   
                                                  999.,//0.98, 
                                                  999.,//1.0, // R9 barrel, R9 endcap   
                                                  0.008, 
                                                  0.008, // Deta barrel, Deta endcap   
                                                  0.07, 
                                                  0.05 // Dphi barrel, Dphi endcap   
                                                  )>=1)   
                {   
                  triggerBit[it] = true;   
                }   
            }   
        }   
    }   
  else if (triggerName.CompareTo("OpenHLT_Ele65_CaloIdVT_TrkIdT")
           == 0)    
    {    
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)    
        {    
          if (prescaleResponse(menu, cfg, rcounter, it))    
            {    
              if (OpenHlt1ElectronSamHarperPassed(65., 0, // ET, L1isolation      
                                                  999., 
                                                  999., // Track iso barrel, Track iso endcap   
                                                  999., 
                                                  999., // Track/pT iso barrel, Track/pT iso endcap   
                                                  999., 
                                                  999., // H/ET iso barrel, H/ET iso endcap   
                                                  999., 
                                                  999., // E/ET iso barrel, E/ET iso endcap   
                                                  0.05, 
                                                  0.05, // H/E barrel, H/E endcap   
                                                  0.011, 
                                                  0.031, // cluster shape barrel, cluster shape endcap   
                                                  0.98, 
                                                  1.0, // R9 barrel, R9 endcap   
                                                  0.008, 
                                                  0.008, // Deta barrel, Deta endcap   
                                                  0.07, 
                                                  0.05 // Dphi barrel, Dphi endcap   
                                                  )>=1)    
                {    
                  triggerBit[it] = true;    
                }    
            }    
        }    
    }    

  // end JH

  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(27., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.125,
						  0.075, // Track/pT iso barrel, Track/pT iso endcap  
						  0.125,
						  0.075, // H/ET iso barrel, H/ET iso endcap  
						  0.125,
						  0.075, // E/ET iso barrel, E/ET iso endcap  
						  0.05,
						  0.05, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(32., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.125,
						  0.075, // Track/pT iso barrel, Track/pT iso endcap  
						  0.125,
						  0.075, // H/ET iso barrel, H/ET iso endcap  
						  0.125,
						  0.075, // E/ET iso barrel, E/ET iso endcap  
						  0.05,
						  0.05, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
		// 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3") == 0 ||
		   triggerName.CompareTo("OpenHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4") == 0    //two names for historical reasons.
		   )
  {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	  {
		  if (prescaleResponse(menu, cfg, rcounter, it))
		  {
			  if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
												  999.,
												  999., // Track iso barrel, Track iso endcap  
												  0.2,
												  0.2, // Track/pT iso barrel, Track/pT iso endcap  
												  0.2,
												  0.2, // H/ET iso barrel, H/ET iso endcap  
												  0.2,
												  0.2, // E/ET iso barrel, E/ET iso endcap  
												  0.15,
												  0.10, // H/E barrel, H/E endcap  
												  0.011,
												  0.031, // cluster shape barrel, cluster shape endcap  
												  999.,
												  999., // R9 barrel, R9 endcap  
												  0.01,
												  0.01, // Deta barrel, Deta endcap  
												  0.15,
												  0.10 // Dphi barrel, Dphi endcap  
												  )>=1)
			  {
				  triggerBit[it] = true;
			  }
		  }
	  }
  }
	
	
	
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele45_CaloIdVT_TrkIdT_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(45., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.05,
						  0.05, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele90_NoSpikeFilter_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(90., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  999.,
						  999., // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  999.,
						  999., // Deta barrel, Deta endcap  
						  999.,
						  999. // Dphi barrel, Dphi endcap  
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoubleEle33_v1")  //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(33., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  999.,
						  999., // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  999.,
						  999., // Deta barrel, Deta endcap  
						  999.,
						  999. // Dphi barrel, Dphi endcap  
						  )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoubleEle33_CaloIdL_v1")  //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(33., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  999.,
						  999., // Deta barrel, Deta endcap  
						  999.,
						  999. // Dphi barrel, Dphi endcap  
						  )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoubleEle8_CaloIdL_TrkIdVL")  //new
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap  
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap  
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap 
						  999.,
						  999., // R9 barrel, R9 endcap  
						  0.01,
						  0.01, // Deta barrel, Deta endcap  
						  0.15,
						  0.10 // Dphi barrel, Dphi endcap  
						  )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-05-11: promoted to v4. Removed R9 cut.
  else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2ElectronsAsymSamHarperPassed(17., 0, // ET, L1isolation 
						       999.,
						       999., // Track iso barrel, Track iso endcap 
						       999.,
						       999., // Track/pT iso barrel, Track/pT iso endcap 
						       0.2,
						       0.2, // H/ET iso barrel, H/ET iso endcap 
						       0.2,
						       0.2, // E/ET iso barrel, E/ET iso endcap
						       0.15,
						       0.10, // H/E barrel, H/E endcap 
						       0.014,
						       0.035, // cluster shape barrel, cluster shape endcap 
						       999.,//0.98,
						       999., // R9 barrel, R9 endcap 
						       999.,
						       999., // Deta barrel, Deta endcap 
						       999.,
						       999., // Dphi barrel, Dphi endcap 
						       8.,
						       0, // ET, L1isolation 
						       999.,
						       999., // Track iso barrel, Track iso endcap 
						       999.,
						       999., // Track/pT iso barrel, Track/pT iso endcap 
						       0.2,
						       0.2, // H/ET iso barrel, H/ET iso endcap 
						       0.2,
						       0.2, // E/ET iso barrel, E/ET iso endcap
						       0.15,
						       0.10, // H/E barrel, H/E endcap 
						       0.014,
						       0.035, // cluster shape barrel, cluster shape endcap 
						       999.,//0.98,
						       999., // R9 barrel, R9 endcap 
						       999.,
						       999., // Deta barrel, Deta endcap 
						       999.,
						       999. // Dphi barrel, Dphi endcap 
						       )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
		// 2011-06-09: added second name to match that in the May 9 train. 
  else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4") == 0 || 
		   triggerName.CompareTo("OpenHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") == 0
		   )
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2ElectronsAsymSamHarperPassed(17., 0, // ET, L1isolation 
						       999.,
						       999., // Track iso barrel, Track iso endcap 
						       0.2,
						       0.2, // Track/pT iso barrel, Track/pT iso endcap  x
						       0.2,
						       0.2, // H/ET iso barrel, H/ET iso endcap  x
						       0.2,
						       0.2, // E/ET iso barrel, E/ET iso endcap x
						       0.15,
						       0.1, // H/E barrel, H/E endcap x
						       0.011,
						       0.031, // cluster shape barrel, cluster shape endcap  x
						       999.,//0.98,
						       999., // R9 barrel, R9 endcap 
						       0.01,
						       0.01, // Deta barrel, Deta endcap x
						       0.15,
						       0.1, // Dphi barrel, Dphi endcap x
						       8.,
						       0, // ET, L1isolation 
						       999.,
						       999., // Track iso barrel, Track iso endcap 
						       0.2,
						       0.2, // Track/pT iso barrel, Track/pT iso endcap 
						       0.2,
						       0.2, // H/ET iso barrel, H/ET iso endcap x
						       0.2,
						       0.2, // E/ET iso barrel, E/ET iso endcap x
						       0.15,//0.1,
						       0.10,//0.075, // H/E barrel, H/E endcap x
						       0.011,
						       0.031, // cluster shape barrel, cluster shape endcap x
						       999.,//0.98,
						       999., // R9 barrel, R9 endcap 
						       0.01,
						       0.01, // Deta barrel, Deta endcap x
						       0.15,
						       0.1 // Dphi barrel, Dphi endcap x
						       )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_DoubleEle8_CaloIdL_TrkIdVL") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   999.,
						   999., // Track/pT iso barrel, Track/pT iso endcap
						   999.,
						   999., // H/ET iso barrel, H/ET iso endcap
						   999.,
						   999., // E/ET iso barrel, E/ET iso endcap
						   0.15,
						   0.1, // H/E barrel, H/E endcap
						   0.014,
						   0.035, // cluster shape barrel, cluster shape endcap
						   999.,//0.98,
						   999.,//1.0, // R9 barrel, R9 endcap
						   0.01,
						   0.01, // Deta barrel, Deta endcap
						   0.15,
						   0.1 // Dphi barrel, Dphi endcap
						   )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
		//2011-06-09
  else if (triggerName.CompareTo("OpenHLT_DoubleEle8_CaloIdT_TrkIdVL_v1") == 0)
  {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	  {
		  if (prescaleResponse(menu, cfg, rcounter, it))
		  {
			  if (OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation
												   999.,
												   999., // Track iso barrel, Track iso endcap
												   999.,
												   999., // Track/pT iso barrel, Track/pT iso endcap
												   999.,
												   999., // H/ET iso barrel, H/ET iso endcap
												   999.,
												   999., // E/ET iso barrel, E/ET iso endcap
												   0.1,
												   0.75, // H/E barrel, H/E endcap x
												   0.011,
												   0.031, // cluster shape barrel, cluster shape endcap x
												   999.,
												   999., // R9 barrel, R9 endcap
												   0.01,
												   0.01, // Deta barrel, Deta endcap x
												   0.15,
												   0.1 // Dphi barrel, Dphi endcap x
												   )>=2)
			  {
				  triggerBit[it] = true;
			  }
		  }
	  }
  }
	
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt3ElectronsSamHarperPassed(10., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   999.,
						   999., // Track/pT iso barrel, Track/pT iso endcap
						   999.,
						   999., // H/ET iso barrel, H/ET iso endcap
						   999.,
						   999., // E/ET iso barrel, E/ET iso endcap
						   0.15,
						   0.1, // H/E barrel, H/E endcap
						   999.,
						   999., // cluster shape barrel, cluster shape endcap
						   999.,//0.98,
						   999.,//1., // R9 barrel, R9 endcap
						   999.,
						   999., // Deta barrel, Deta endcap
						   999.,
						   999. // Dphi barrel, Dphi endcap
						   )>=3 && OpenHlt2ElectronsSamHarperPassed(10., 0, // ET, L1isolation
											    999.,
											    999., // Track iso barrel, Track iso endcap
											    999.,
											    999., // Track/pT iso barrel, Track/pT iso endcap
											    999.,
											    999., // H/ET iso barrel, H/ET iso endcap
											    999.,
											    999., // E/ET iso barrel, E/ET iso endcap
											    0.15,
											    0.10, // H/E barrel, H/E endcap 
											    0.014,
											    0.035, // cluster shape barrel, cluster shape endcap 
											    0.98,
											    1., // R9 barrel, R9 endcap
											    0.01,
											    0.01, // Deta barrel, Deta endcap 
											    0.15,
											    0.10 // Dphi barrel, Dphi endcap
											    )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_TripleEle10_CaloIdL_TrkIdVL_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt3ElectronsSamHarperPassed(10., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   999.,
						   999., // Track/pT iso barrel, Track/pT iso endcap
						   999.,
						   999., // H/ET iso barrel, H/ET iso endcap
						   999.,
						   999., // E/ET iso barrel, E/ET iso endcap
						   0.15,
						   0.10, // H/E barrel, H/E endcap 
						   0.014,
						   0.035, // cluster shape barrel, cluster shape endcap 
						   999.,//0.98,
						   999.,//1., // R9 barrel, R9 endcap 
						   0.01,
						   0.01, // Deta barrel, Deta endcap 
						   0.15,
						   0.10 // Dphi barrel, Dphi endcap
						   )>=3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // proxy to the HLT implementation TODO should use SC variable once available in ntuples
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(
						  32., 
						  0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999,
						  999, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.10, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  999,
						  999, // Deta barrel, Deta endcap  
						  999,
						  999 // Dphi barrel, Dphi endcap  
						  ) >=1
		  && 
		  OpenHlt1PhotonSamHarperPassed(
						17., 
						0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.1, // H/E barrel, H/E endcap
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999,
						999, // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						) >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(
						  32.,                         
						  0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  0.125,
						  0.075, // Track/pT iso barrel, Track/pT iso endcap  
						  0.125,
						  0.075, // H/ET iso barrel, H/ET iso endcap  
						  0.125,
						  0.075, // E/ET iso barrel, E/ET iso endcap  
						  0.10,
						  0.075, // H/E barrel, H/E endcap  
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap  
						  999.,
						  999., // R9 barrel, R9 endcap  
						  0.008,
						  0.008, // Deta barrel, Deta endcap  
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap  
						  ) >=1
		  && 
		  OpenHlt1PhotonSamHarperPassed(
						17., 
						0, // ET, L1isolation  
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.1, // H/E barrel, H/E endcap 
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999,
						999, // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						) >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-05-11: promoted to v4, removed R9, and edited mass filtration section to include identical Et exclusion 
  else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1ElectronSamHarperPassed(
										   17., 
										   0, // ET, L1isolation  
										   999.,
										   999., // Track iso barrel, Track iso endcap  
										   0.05,
										   0.05, // Track/pT iso barrel, Track/pT iso endcap  
										   0.05,
										   0.05, // H/ET iso barrel, H/ET iso endcap  
										   0.05,
										   0.05, // E/ET iso barrel, E/ET iso endcap  
										   0.05,
										   0.05, // H/E barrel, H/E endcap  
										   0.011,
										   0.031, // cluster shape barrel, cluster shape endcap  
										   999.,//0.98,
										   999.,//1.0, // R9 barrel, R9 endcap  
										   0.008,
										   0.008, // Deta barrel, Deta endcap  
										   0.07,
										   0.05 );// Dphi barrel, Dphi endcap 
	      if (firstVector.size()>=1){
		std::vector<int> secondVector = VectorOpenHlt1PhotonSamHarperPassed(
										    8., 
										    0, // ET, L1isolation
										    999.,
										    999., // Track iso barrel, Track iso endcap
										    999.,
										    999., // Track/pT iso barrel, Track/pT iso endcap
										    999.,
										    999., // H iso barrel, H iso endcap
										    999.,
										    999., // E iso barrel, E iso endcap
										    0.15,
										    0.1, // H/E barrel, H/E endcap
										    999.,
										    999., // cluster shape barrel, cluster shape endcap
										    999.,
										    999., // R9 barrel, R9 endcap
										    999.,
										    999., // Deta barrel, Deta endcap
										    999.,
										    999.); // Dphi barrel, Dphi endcap
		if (secondVector.size()>=2){
						
						
		  // mass condition
		  TLorentzVector ele;
		  TLorentzVector pho;
		  TLorentzVector sum;
		  float mass = 0.;
		  for (unsigned int i=0; i<firstVector.size(); i++)
		    {
                      ele.SetPtEtaPhiM(
                                       ohEleEt[firstVector[i]],
                                       ohEleEta[firstVector[i]],
                                       ohElePhi[firstVector[i]],
                                       0.);

		      for (unsigned int j=0; j<secondVector.size() ; j++)
			{
								
			  if(ohEleEt[firstVector[i]] - ohPhotEt[secondVector[j]] < 0.00001) continue;
			  pho.SetPtEtaPhiM(
					   ohPhotEt[secondVector[j]],
					   ohPhotEta[secondVector[j]],
					   ohPhotPhi[secondVector[j]],
					   0.);
			  sum = ele + pho;
			  mass = sum.M();
								
			  if (mass>=30)
			    triggerBit[it] = true;
								
			}//endfor
		    }// end for
						
		}//end if size>=2
	      }//end if first vector size >= 1
	    }
	}
    }//end trigger
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1ElectronSamHarperPassed(
										   17., //
										   0, // ET, L1isolation  
										   999.,
										   999., // Track iso barrel, Track iso endcap  
										   0.05,
										   0.05, // Track/pT iso barrel, Track/pT iso endcap  
										   0.05,
										   0.05, // H/ET iso barrel, H/ET iso endcap  
										   0.05,
										   0.05, // E/ET iso barrel, E/ET iso endcap  
										   0.05,
										   0.05, // H/E barrel, H/E endcap  
										   0.011,
										   0.031, // cluster shape barrel, cluster shape endcap  
										   999.,
										   999., // R9 barrel, R9 endcap  //this was broke in full hlt
										   0.008,
										   0.008, // Deta barrel, Deta endcap  
										   0.07,
										   0.05 );// Dphi barrel, Dphi endcap  
	      if (firstVector.size()>=1){
		std::vector<int> secondVector = VectorOpenHlt1ElectronSamHarperPassed(
										      8., 
										      0, // ET, L1isolation
										      999.,
										      999., // Track iso barrel, Track iso endcap
										      999.,
										      999., // Track/pT iso barrel, Track/pT iso endcap
										      999.,
										      999., // H iso barrel, H iso endcap
										      999.,
										      999., // E iso barrel, E iso endcap
										      0.15,
										      0.1, // H/E barrel, H/E endcap 
										      999.,
										      999., // cluster shape barrel, cluster shape endcap
										      999.,
										      999., // R9 barrel, R9 endcap
										      999.,
										      999., // Deta barrel, Deta endcap
										      999.,
										      999.); // Dphi barrel, Dphi endcap
		if (secondVector.size()>=2){
						
						
		  // mass condition
		  TLorentzVector ele;
		  TLorentzVector pho;
		  TLorentzVector sum;
		  float mass = 0.;
		  for (unsigned int i=0; i<firstVector.size(); i++)
		    {
		      for (unsigned int j=0; j<secondVector.size() ; j++)
			{
								
			  //                  if (firstVector[i] == secondVector[j]) continue;
			  ele.SetPtEtaPhiM(
					   ohEleEt[firstVector[i]],
					   ohEleEta[firstVector[i]],
					   ohElePhi[firstVector[i]],
					   0.);
			  pho.SetPtEtaPhiM(
					   ohPhotEt[secondVector[j]],
					   ohPhotEta[secondVector[j]],
					   ohPhotPhi[secondVector[j]],
					   0.);
			  sum = ele + pho;
			  mass = sum.M();
								
			  if (mass>30)
			    triggerBit[it] = true;
								
			}//endfor
		    }// end for
						
		}//end if size>=2
	      }//end if first vector size >= 1
	    }
	}
    }//end trigger
	
  /* Photons */
	
else if (triggerName.CompareTo("OpenHLT_Photon30_CaloIdVT_CentralJet20_BTagIP") == 0 )
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            //return the highest Pt photon id, return -999, if there is no qualified photon
	   std::vector<int> photonIndicesVector = VectorOpenHlt1PhotonSamHarperPassed(30., 0, // ET, L1isolation
                  999,
                  999, // Track iso barrel, Track iso endcap
                  999,
                  999, // Track/pT iso barrel, Track/pT iso endcap
                  999,
                  999, // H iso barrel, H iso endcap
                  999,
                  999, // E iso barrel, E iso endcap
                  0.05,
                  0.05, // H/E barrel, H/E endcap
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap
                  999,
                  999, // R9 barrel, R9 endcap
                  999.,
                  999., // Deta barrel, Deta endcap
                  999.,
                  999. // Dphi barrel, Dphi endcap
										      ) ;
	     if (photonIndicesVector.size()>=1){

          
	       for (int j = 0; j < NohBJetL2Corrected; j++)
		 {//loop over jets
		   if (ohBJetL2CorrectedEt[j] > 20. && fabs(ohBJetL2CorrectedEta[j]) < 2.6)
		     {// ET and eta cuts
		       //                        if (ohBJetIPL25Tag[j] >= 0)
		     if (true)
		       {// Level 2.5 b tag
			 if (ohBJetIPL3Tag[j] >= 3.3)
			   {// Level 3 b tag
			     double deltaEta = fabs(ohBJetL2CorrectedEta[j]-ohPhotEta[photonIndicesVector[0]]);
			     double deltaPhi = fabs(ohBJetL2CorrectedPhi[j]-ohPhotPhi[photonIndicesVector[0]]) < 3.14159 ? ( fabs(ohBJetL2CorrectedPhi[j]-ohPhotPhi[photonIndicesVector[0]])):( fabs(ohBJetL2CorrectedPhi[j]-ohPhotPhi[photonIndicesVector[0]]) - 3.14159);
			     double deltaR = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
			     if (deltaR > 0)
			       {
				 triggerBit[it] = true;
				 break;
			       }
			   }
		       }
		   }// ET and eta cuts
		 }//end loop over jets
	     }
	 }
      }
   }


  else if (triggerName.CompareTo("OpenHLT_Photon20_CaloIdVL_IsoL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(20., 0, // ET, L1isolation
						3.5,
						3.5, // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						3.5,
						3.5, // H iso barrel, H iso endcap
						5.5,
						5.5, // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap
						999.,//0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-05-11 promoted to v4. Removed R9 cut.
  else if (triggerName.CompareTo("OpenHLT_Photon30_CaloIdVL_v4") == 0){
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	  {
		  if (prescaleResponse(menu, cfg, rcounter, it))
		  {
			  if (OpenHlt1PhotonSamHarperPassed(30., 0, // ET, L1isolation
												999.,
												999., // Track iso barrel, Track iso endcap
												999.,
												999., // Track/pT iso barrel, Track/pT iso endcap
												999.,
												999., // H iso barrel, H iso endcap
												999.,
												999., // E iso barrel, E iso endcap
												0.15,
												0.10, // H/E barrel, H/E endcap 
												0.024,
												0.040, // cluster shape barrel, cluster shape endcap 
												999.,//0.98,
												999., // R9 barrel, R9 endcap
												999.,
												999., // Deta barrel, Deta endcap
												999.,
												999. // Dphi barrel, Dphi endcap
												)>=1)
			  {
				  triggerBit[it] = true;
			  }
		  }
	  }
  }
/*
//2011-6-8
  else if (triggerName.CompareTo("OpenHLT_Photon10_CaloIdVL_v3") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(10., 0, // ET, L1isolation
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                0.15,
                                                0.10, // H/E barrel, H/E endcap 
                                                0.024,
                                                0.040, // cluster shape barrel, cluster shape endcap 
                                                999.,//0.98,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }

//2011-6-8
  else if (triggerName.CompareTo("OpenHLT_Photon15_CaloIdVL_v3") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(15., 0, // ET, L1isolation
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                0.15,
                                                0.10, // H/E barrel, H/E endcap 
                                                0.024,
                                                0.040, // cluster shape barrel, cluster shape endcap 
                                                999.,//0.98,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
*/

//end new	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon50_CaloIdVL_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(50., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap  
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap   
						999.,
						999., // R9 barrel, R9 endcap  
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-05-11 promoted to v4 
  else if (triggerName.CompareTo("OpenHLT_Photon30_CaloIdVL_IsoL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(30., 0, // ET, L1isolation
						3.5,
						3.5, // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						3.5,
						3.5, // H iso barrel, H iso endcap
						5.5,
						5.5, // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap
						999.,//0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	//2011-05-11 promoted to v3, removed R9 cut.
  else if (triggerName.CompareTo("OpenHLT_Photon50_CaloIdVL_IsoL_v3") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(50., 0, // ET, L1isolation
						3.5,
						3.5, // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						3.5,
						3.5, // H iso barrel, H iso endcap
						5.5,
						5.5, // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap
						999.,//0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-05-11 promoted to v4, removed R9 filter. 
  else if (triggerName.CompareTo("OpenHLT_Photon75_CaloIdVL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(75., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap
						999.,//0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon90_CaloIdVL_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(90., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap 
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap 
						999.,
						999., // R9 barrel, R9 endcap 
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-05-11 promoted to v1, removed R9 cut. 
  else if (triggerName.CompareTo("OpenHLT_Photon75_CaloIdVL_IsoL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(75., 0, // ET, L1isolation
						3.5,
						3.5, // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						3.5,
						3.5, // H iso barrel, H iso endcap
						5.5,
						5.5, // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap
						999.,//0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon90_CaloIdVL_IsoL_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(90., 0, // ET, L1isolation
						3.5,
						3.5, // Track iso barrel, Track iso endcap 
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						3.5,
						3.5, // H iso barrel, H iso endcap 
						5.5,
						5.5, // E iso barrel, E iso endcap 
						0.15,
						0.10, // H/E barrel, H/E endcap 
						0.024,
						0.040, // cluster shape barrel, cluster shape endcap 
						999.,
						999., // R9 barrel, R9 endcap 
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Photon125_NoSpikeFilter_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(125., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						999.,
						999., // H/E barrel, H/E endcap
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon125_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(125., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap 
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap 
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon200_NoHE_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(200., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						999.,
						999., // H/E barrel, H/E endcap
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }	
  //2011-06-08
  else if (triggerName.CompareTo("OpenHLT_Photon225_NoHE_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(225., 0, // ET, L1isolation 
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                999.,
                                                999., // H/E barrel, H/E endcap
                                                999.,
                                                999., // cluster shape barrel, cluster shape endcap
                                                999.,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }  

	
  // 2011-05-11 promoted to v4 
  else if (triggerName.CompareTo("OpenHLT_Photon32_CaloIdL_Photon26_CaloIdL_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
                //That's OpenHLTPhoCuts(e_t, H/E_EB, H/E_EC, cluster_shape_EB, cluster_shape_EC, cal_iso/Et, trk_iso/Et)
                //cal_iso is both HCAL_iso and ECAL_iso
	      if (OpenHltPhoCuts(32, 0.15, 0.10, 0.014, 0.035, 999., 999.) >= 1
		  && OpenHltPhoCuts(26, 0.15, 0.10, 0.014, 0.035, 999., 999.)
		  >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-05-11 promoted to v3
  else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdL_Photon22_CaloIdL_v3") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
		//That's OpenHLTPhoCuts(e_t, H/E_EB, H/E_EC, cluster_shape_EB, cluster_shape_EC, cal_iso/Et, trk_iso/Et)
		//cal_iso is both HCAL_iso and ECAL_iso
	      if (OpenHltPhoCuts(36, 0.15, 0.10, 0.014, 0.035, 999., 999.) >= 1
		  && OpenHltPhoCuts(22, 0.15, 0.10, 0.014, 0.035, 999., 999.)
		  >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

//2011-05-22
   else if (triggerName.CompareTo("OpenHLT_Photon48_CaloIdL_Photon38_CaloIdL_v1") == 0)
   {
	   if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	   {
		   if (prescaleResponse(menu, cfg, rcounter, it))
		   {
			   if (OpenHltPhoCuts(48, 0.15, 0.10, 0.014, 0.035, 999, 999) >= 1
				   && OpenHltPhoCuts(38, 0.15, 0.10, 0.014, 0.035, 999, 999)
				   >= 2)
			   {
				   triggerBit[it] = true;
			   }
		   }
	   }
   }
	
	//2011-05-22
   else if (triggerName.CompareTo("OpenHLT_Photon40_CaloIdL_Photon28_CaloIdL_v1") == 0)
   {
	   if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	   {
		   if (prescaleResponse(menu, cfg, rcounter, it))
		   {
			   if (OpenHltPhoCuts(40, 0.15, 0.10, 0.014, 0.035, 999, 999) >= 1
				   && OpenHltPhoCuts(28, 0.15, 0.10, 0.014, 0.035, 999, 999)
				   >= 2)
			   {
				   triggerBit[it] = true;
			   }
		   }
	   }
   }

		//2011-06-09
   else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdVL_Photon22_CaloIdVL_v1") == 0)
   {
	   if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	   {
		   if (prescaleResponse(menu, cfg, rcounter, it))
		   {
			   if (OpenHltPhoCuts(36, 0.15, 0.10, 0.024, 0.040, 999, 999) >= 1
				   && OpenHltPhoCuts(22, //et
									 0.15, // h/e_eb
									 0.10, // h/e_ec 
									 0.024, //cluster shape eb
									 0.040, //cluster shape ec
									 999, //cal_iso
									 999) //tra_iso
				   >= 2)
			   {
				   triggerBit[it] = true;
			   }
		   }
	   }
   }
	
	
//2011-06-08
   else if (triggerName.CompareTo("OpenHLT_Photon44_CaloIdL_Photon34_CaloIdL_v1") == 0)
   {
           if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
           {
                   if (prescaleResponse(menu, cfg, rcounter, it))
                   {
                           if (OpenHltPhoCuts(44, 0.15, 0.10, 0.014, 0.035, 999, 999) >= 1
							   && OpenHltPhoCuts(34, //et
												 0.15, // h/e_eb
												 0.10, // h/e_ec 
												 0.014, //cluster shape eb
												 0.035, //cluster shape ec
												 999, //cal_iso
												 999) //tra_iso
							   >= 2)
                           {
                                   triggerBit[it] = true;
                           }
                   }
           }
   }

	
  // to be removed ??
  /*PhotonX_(M)HTX */
	
  else if (isPhotonX_HTXTrigger(triggerName, thresholds, r9Id, caloId, photonIso))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonPassed(thresholds[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]
				     ) >= 1	&& OpenHltSumCorHTPassed(thresholds[1], 40.)>=1)
	      
	      {
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
   else if (isPhotonX_MHTXTrigger(triggerName, thresholds, r9Id, caloId, photonIso))
     {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonPassed(thresholds[0],
				     map_PhotonR9ID[r9Id[0]],
				     map_EGammaCaloId[caloId[0]],
				     map_PhotonIso[photonIso[0]]
				     ) >= 1  && OpenHltMHT(thresholds[1], 30.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
 	}
     }
	
//   else if (isPhotonX_CaloIdL_HTXTrigger(triggerName, thresholds))
//     {
//       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
// 	{
// 	  if (prescaleResponse(menu, cfg, rcounter, it))
// 	    {
// 	      if (OpenHltPhoCuts(thresholds[0], 0.15, 0.10, 0.014, 0.034, 999, 999) >= 1
// 		  && OpenHltSumCorHTPassed(thresholds[1], 30.)>=1)
// 		{
// 		  triggerBit[it] = true;
// 		}
// 	    }
// 	}
//     }
	
//   else if (isPhotonX_CaloIdL_MHTXTrigger(triggerName, thresholds))
//     {
//       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
// 	{
// 	  if (prescaleResponse(menu, cfg, rcounter, it))
// 	    {
// 	      if (OpenHltPhoCuts(thresholds[0], 0.15, 0.10, 0.014, 0.034, 999, 999) >= 1
// 		  && OpenHltMHT(thresholds[1], 30.)>=1)
// 		{
// 		  triggerBit[it] = true;
// 		}
// 	    }
// 	}
//     }
	
  /*
   * 2011-01-26 added by Christian Hartl
   * according to https://twiki.cern.ch/twiki/bin/view/CMS/EgammaWorkingPoints
   * adapted from "OpenHLT_Photon65_CaloEleId_Isol_L1R":
   * - dropped _L1R from name because L1 not relaxed anymore 
   * TODO
   * - Tisobarrel=0.001, Tisoendcap=5.0 ? -- check with Mass.
   */
  // to be removed
//   else if (triggerName.CompareTo("OpenHLT_Photon65_CaloEleId_Isol") == 0)
//     {
//       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
// 	{
// 	  if (prescaleResponse(menu, cfg, rcounter, it))
// 	    {
// 	      if (OpenHlt1PhotonSamHarperPassed( 65., 0, // ET, L1isolation
// 						 0.001,
// 						 5.0, // Track iso barrel, Track iso endcap ???
// 						 0.2,
// 						 0.1, // Track/pT iso barrel, Track/pT iso endcap
// 						 0.2,
// 						 0.1, // H iso barrel, H iso endcap
// 						 0.2,
// 						 0.1, // E iso barrel, E iso endcap
// 						 0.05,
// 						 0.05, // H/E barrel, H/E endcap
// 						 0.014,
// 						 0.035, // cluster shape barrel, cluster shape endcap
// 						 999.,//0.98,
// 						 999.,
// 						 999., // Deta barrel, Deta endcap
// 						 999.,
// 						 999. // Dphi barrel, Dphi endcap
// 						 )>=1)
// 		{
// 		  triggerBit[it] = true;
// 		}
// 	    }
// 	}
//     }
	
	
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton32_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(32., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						0.014,
						0.035, // cluster shape barrel, cluster shape endcap
						0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton33_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(33., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap
						999.,
						999., // cluster shape barrel, cluster shape endcap
						0.98,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton33_HEVT_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(33., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.05,
						0.05, // H/E barrel, H/E endcap   
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap 
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011 June 08, Arnaud , copied from above
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton38_HEVT_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(38., 0, // ET, L1isolation
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.05,
						0.05, // H/E barrel, H/E endcap   
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap 
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton50_v1") == 0) //new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(50., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap 
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_DoublePhoton60_v1") == 0) //new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PhotonSamHarperPassed(60., 0, // ET, L1isolation 
						999.,
						999., // Track iso barrel, Track iso endcap
						999.,
						999., // Track/pT iso barrel, Track/pT iso endcap
						999.,
						999., // H iso barrel, H iso endcap
						999.,
						999., // E iso barrel, E iso endcap
						0.15,
						0.10, // H/E barrel, H/E endcap 
						999.,
						999., // cluster shape barrel, cluster shape endcap
						999.,
						999., // R9 barrel, R9 endcap
						999.,
						999., // Deta barrel, Deta endcap
						999.,
						999. // Dphi barrel, Dphi endcap
						)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

// 2011-06-08
 else if (triggerName.CompareTo("OpenHLT_DoublePhoton80_v1") == 0) //new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(80., 0, // ET, L1isolation 
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                0.15,
                                                0.10, // H/E barrel, H/E endcap 
                                                999.,
                                                999., // cluster shape barrel, cluster shape endcap
                                                999.,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=2)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
	
//new
// 2011-06-08
 else if (triggerName.CompareTo("OpenHLT_Photon135_v1") == 0) //new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(135., 0, // ET, L1isolation 
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                0.15,
                                                0.10, // H/E barrel, H/E endcap 
                                                999.,
                                                999., // cluster shape barrel, cluster shape endcap
                                                999.,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }


// 2011-06-08
 else if (triggerName.CompareTo("OpenHLT_Photon400_v1") == 0) //new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1PhotonSamHarperPassed(400., 0, // ET, L1isolation 
                                                999.,
                                                999., // Track iso barrel, Track iso endcap
                                                999.,
                                                999., // Track/pT iso barrel, Track/pT iso endcap
                                                999.,
                                                999., // H iso barrel, H iso endcap
                                                999.,
                                                999., // E iso barrel, E iso endcap
                                                0.15,
                                                0.10, // H/E barrel, H/E endcap 
                                                999.,
                                                999., // cluster shape barrel, cluster shape endcap
                                                999.,
                                                999., // R9 barrel, R9 endcap
                                                999.,
                                                999., // Deta barrel, Deta endcap
                                                999.,
                                                999. // Dphi barrel, Dphi endcap
                                                )>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }


	
  /* Taus */
  else if (triggerName.CompareTo("OpenHLT_DoubleIsoTau15_Trk5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltTauL2SCPassed(15., 5., 0, 0., 1, 14., 30.)>=2)
		{ //Thresholds are for UNcorrected L1 jets in 8E29
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isIsoPFTauX_TrkX_METXTrigger(triggerName, thresholds)){
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1){
      if (prescaleResponse(menu, cfg, rcounter, it)){
        if (OpenHltTightConeIsoPFTauPassed(thresholds[0], 99.,  thresholds[1], 14., 30.)>=1 ){
	  if (recoMetCal > thresholds[2]) {
	    triggerBit[it] = true;
	  }
        }
      }
    }
  }

  else if (isIsoPFTauX_TrkXTrigger(triggerName, thresholds)){ 
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1){ 
      if (prescaleResponse(menu, cfg, rcounter, it)){ 
	if (OpenHltTightConeIsoPFTauPassed(thresholds[0], 99.,  thresholds[1], 14., 30.)>=1 ){ 
	  triggerBit[it] = true; 
	} 
      }
    } 
  }

  else if (isDoubleIsoPFTauX_X_TrkX_eta2pXTrigger(triggerName, thresholds)){
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1){
      if (prescaleResponse(menu, cfg, rcounter, it)){
	if (OpenHltTightConeIsoPFTauPassed(thresholds[0], thresholds[3],  thresholds[2], 14., 30.)>=1 &&
	    OpenHltTightConeIsoPFTauPassed(thresholds[1], thresholds[3],  thresholds[2], 14., 30.)>=2 ){ 
	  triggerBit[it] = true;
	}
      }
    }
  }

  else if (isDoubleLooseIsoPFTauX_X_TrkX_eta2pXTrigger(triggerName, thresholds)){
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1){
      if (prescaleResponse(menu, cfg, rcounter, it)){
	if (OpenHltIsoPFTauPassed(thresholds[0], thresholds[3],  thresholds[2], 14., 30.)>=1 &&
	    OpenHltIsoPFTauPassed(thresholds[1], thresholds[3],  thresholds[2], 14., 30.)>=2 ){ 
	  triggerBit[it] = true;
	}
      }
    }
  }

  /* BTag */
	
  /**********BTagMu_JetX(U)**********/
  else if (isBTagMu_JetXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
				
	      // apply L2 cut on jets 
	      for (int i = 0; i < NohBJetL2; i++)
		if (ohBJetL2Et[i] > thresholds[0] && fabs(ohBJetL2Eta[i]) < 3.0) // change this ET cut to 20 for the 20U patath 
		  njets++;
				
	      // apply b-tag cut 
	      int max = (NohBJetL2 > 4) ? 4 : NohBJetL2;
	      for (int i = 0; i < max; i++)
		{
		  if (ohBJetL2Et[i] > 10.)
		    { // keep this at 10 even for the 20UU path - also, no eta cut here 
		      if (ohBJetPerfL25Tag[i] > 0.5)
			{ // Level 2.5 b tag 
			  if (ohBJetPerfL3Tag[i] > 0.5)
			    { // Level 3 b tag 
			      rc++;
			    }
			}
		    }
		}
	      if (rc >= 1 && njets>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (isBTagMu_JetXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int max = (NohBJetL2Corrected > 2) ? 2 : NohBJetL2Corrected;
	      for (int i = 0; i < max; i++)
		{
		  if (ohBJetL2CorrectedEt[i] > thresholds[0])
		    { // ET cut
		      if (ohBJetPerfL25Tag[i] > 0.5)
			{ // Level 2.5 b tag
			  if (ohBJetPerfL3Tag[i] > 0.5)
			    { // Level 3 b tag
			      rc++;
			    }
			}
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /******BTagMu_DiJetXU******/
	
  else if (isBTagMu_DiJetXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
				
	      // apply L2 cut on jets
	      for (int i = 0; i < NohBJetL2; i++)
		if (ohBJetL2Et[i] > thresholds[0] && fabs(ohBJetL2Eta[i]) < 3.0) // change this ET cut to 20 for the 20U patath
		  njets++;
				
	      // apply b-tag cut
	      for (int i = 0; i < NohBJetL2; i++)
		{
		  if (ohBJetL2Et[i] > 10.)
		    { // keep this at 10 even for the 20UU path - also, no eta cut here
		      if (ohBJetPerfL25Tag[i] > 0.5)
			{ // Level 2.5 b tag
			  if (ohBJetPerfL3Tag[i] > 0.5)
			    { // Level 3 b tag
			      rc++;
			    }
			}
		    }
		}
	      if (rc >= 1 && njets>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  /**********BTagMu_DiJetXU_MuX***************************************/
  else if (isBTagMu_DiJetXU_MuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
				
	      // apply L2 cut on jets
	      for (int i = 0; i < NohBJetL2; i++)
		if (ohBJetL2Et[i] > thresholds[0] && fabs(ohBJetL2Eta[i]) < 3.0)
		  njets++;
				
	      // apply b-tag cut
	      for (int i = 0; i < NohBJetL2; i++)
		{
		  if (ohBJetL2Et[i] > 10.)
		    { // keep this at 10 even for all btag mu paths
		      if (ohBJetPerfL25Tag[i] > 0.5)
			{ // Level 2.5 b tag
			  if (OpenHlt1L3MuonPassed(thresholds[1], 5.0) >=1)
			    {//require at least one L3 muon
			      if (ohBJetPerfL3Tag[i] > 0.5)
				{ // Level 3 b tag
				  rc++;
				}
			    }
			}
		    }
		}
	      if (rc >= 1 && njets>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /**********Available in 2011 menu: BTagMu_DiJetX_MuX***************************************/
  else if (isBTagMu_DiJetX_MuXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
				
	      // apply L2 cut on jets
	      for (int i = 0; i < NohBJetL2Corrected; i++)
		if (ohBJetL2CorrectedEt[i] > thresholds[0] && fabs(ohBJetL2Eta[i]) < 3.0)
		  njets++;
				
	      // apply b-tag cut
	      for (int i = 0; i < NohBJetL2Corrected; i++)
		{
		  if (ohBJetL2CorrectedEt[i] > thresholds[0])
		    { // keep this at 10 even for all btag mu paths
		      if (ohBJetPerfL25Tag[i] > 0.5)
			{ // Level 2.5 b tag
			  if (OpenHlt1L3MuonPassed(thresholds[1], 5.0) >=1)
			    {//require at least one L3 muon
			      if (ohBJetPerfL3Tag[i] > 0.5)
				{ // Level 3 b tag
				  rc++;
				}
			    }
			}
		    }
		}
	      if (rc >= 1 && njets>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /****************BTagIP_JetX*********************************/
	
  else if (isBTagIP_JetXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int max = (NohBJetL2Corrected > 2) ? 2 : NohBJetL2Corrected;
	      for (int i = 0; i < max; i++)
		{
		  if (ohBJetL2CorrectedEt[i] > thresholds[0])
		    { // ET cut 
		      if (ohBJetIPL25Tag[i] > 2.5)
			{ // Level 2.5 b tag 
			  if (ohBJetIPL3Tag[i] > 3.5)
			    { // Level 3 b tag 
			      rc++;
			    }
			}
		    }
		}
	      if (rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /****************QuadJetX_BTagIP*********************************/
	
  else if (isQuadJetX_BTagIPTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltQuadCorJetPassed(thresholds[0])>=1)
		{
		  int rc = 0;
		  int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
		  for (int i = 0; i < max; i++)
		    {
		      if (ohBJetL2CorrectedEt[i] > thresholds[0])
			{ // ET cut 
			  if (ohBJetIPL25Tag[i] > 0)
			    { // Level 2.5 b tag 
			      if (ohBJetIPL3Tag[i] > 2.0)
				{ // Level 3 b tag 
				  rc++;
				}
			    }
			}
		    }
		  if (rc >= 1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
  //Hbb mu+jet trigger
  else if (triggerName.CompareTo("OpenHLT_Mu12_DiCentralJet30_BTagIP3D") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++)
		if (ohBJetL2CorrectedEt[j] > 30. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  njets ++;
	      
	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 30.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 2.5)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 3.5)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && rc >= 1 && OpenHlt1MuonPassed(0., 10., 12., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

 else if (triggerName.CompareTo("OpenHLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      int njetsht = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++){
		if (ohBJetL2CorrectedEt[j] > 38. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  {
		    njets ++;
		    if (ohBJetL2CorrectedEt[j] > 46.) njetsht++;
		  }
	      }

	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 4.0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 6.0)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && njetsht >= 1 && rc >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  //renamed trigger
 else if (triggerName.CompareTo("OpenHLT_CentralJet46_CentralJet38_DiBTagIP3D") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      int njetsht = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++){
		if (ohBJetL2CorrectedEt[j] > 38. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  {
		    njets ++;
		    if (ohBJetL2CorrectedEt[j] > 46.) njetsht++;
		  }
	      }

	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 3.0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 6.0)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && njetsht >= 1 && rc >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  //high pT trigger
 else if (triggerName.CompareTo("OpenHLT_CentralJet60_CentralJet53_DiBTagIP3D") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      int njetsht = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++){
		if (ohBJetL2CorrectedEt[j] > 60. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  {
		    njets ++;
		    if (ohBJetL2CorrectedEt[j] > 53.) njetsht++;
		  }
	      }

	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 2.5)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 3.5)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && njetsht >= 1 && rc >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  //July7th bPhi triggers
 else if (triggerName.CompareTo("OpenHLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      int njetsht = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++){
		if (ohBJetL2CorrectedEt[j] > 38. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  {
		    njets ++;
		    if (ohBJetL2CorrectedEt[j] > 46.) njetsht++;
		  }
	      }

	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 3.0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 6.0)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && njetsht >= 1 && rc >= 2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  //control region trigger
 else if (triggerName.CompareTo("OpenHLT_DiCentralJet36_BTagIP3DLoose") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
	      int njets = 0;
	      
	      for (int j = 0; j < NohBJetL2Corrected ; j++){
		if (ohBJetL2CorrectedEt[j] > 36. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts
		  {
		    njets ++;
		  }
	      }

	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
		  
		  if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets
		    
		    if (ohBJetIPL25Tag[j] > 0.)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] > 3.0)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET cut
		}//loop over jets
				
	      if (njets >=2 && rc >= 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

 else if (triggerName.CompareTo("OpenHLT_DiCentralJet46_BTagIP3DLoose") == 0) 
   { 
     if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
       { 
	 if (prescaleResponse(menu, cfg, rcounter, it)) 
	   { 
	     int rc = 0; 
	     int njets = 0; 
               
	     for (int j = 0; j < NohBJetL2Corrected ; j++){ 
	       if (ohBJetL2CorrectedEt[j] > 46. && fabs(ohBJetL2CorrectedEta[j]) < 2.6) // ET and eta cuts 
		 { 
		   njets ++; 
		 } 
	     } 
 
	     int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected; 
	     for (int j = 0; j < max; j++) 
	       {//loop over jets 
                   
		 if (ohBJetL2CorrectedEt[j] > 20.) {// ET cut on L2.5 jets 
                     
		   if (ohBJetIPL25Tag[j] > 0.) 
		     { // Level 2.5 b tag   
		       if (ohBJetIPL3Tag[j] > 3.0) 
			 { // Level 3 b tag   
			   rc++; 
			 } 
		     } 
		 } // ET cut 
	       }//loop over jets 
                                 
	     if (njets >=2 && rc >= 1) 
	       { 
		 triggerBit[it] = true; 
	       } 
	   } 
       } 
   } 


  /**********************************************/
	
  // Single Top Triggers : Muon channel
	
  /**********************************************/
		
  else if (triggerName.CompareTo("OpenHLT_Mu17_CentralJet30_BTagIP") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
				
	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
					
		  if (ohBJetL2CorrectedEt[j] > 30. && fabs(ohBJetL2CorrectedEta[j]) < 3.) {// ET and eta cuts
						
		    if (ohBJetIPL25Tag[j] >= 0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] >= 3.3)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET and eta cuts
		}//loop over jets
				
	      if (rc >= 1 && OpenHlt1MuonPassed(0., 10., 17., 2., 0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /**********************************************/
	
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_CentralJet30_BTagIP") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
				
	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)	     
		{//loop over jets
					
		  if (ohBJetL2CorrectedEt[j] > 30. && fabs(ohBJetL2CorrectedEta[j]) < 3.) { // ET and eta cuts
						
		    if (ohBJetIPL25Tag[j] >= 0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] >= 3.3)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET and eta cuts
		}//loop over jets
				
	      if (rc >= 1 && OpenHlt1MuonPassed(0., 10., 17., 2., 1)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  /**********************************************/
  // with L2Q and eta restrictions
  else if (triggerName.CompareTo("OpenHLT_Mu17_eta2p1_CentralJet30_BTagIP") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
				
	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)
		{//loop over jets
					
		  if (ohBJetL2CorrectedEt[j] > 30. && fabs(ohBJetL2CorrectedEta[j]) < 3.) {// ET and eta cuts
						
		    if (ohBJetIPL25Tag[j] >= 0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] >= 3.3)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET and eta cuts
		}//loop over jets
				
	      if (rc >= 1 && OpenHlt1MuonPassed(0., 14., 17., 2., 0, 2.1, 2.1, 1, 2)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /**********************************************/
	
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1_CentralJet30_BTagIP") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      int rc = 0;
				
	      int max = (NohBJetL2Corrected > 4) ? 4 : NohBJetL2Corrected;
	      for (int j = 0; j < max; j++)	     
		{//loop over jets
					
		  if (ohBJetL2CorrectedEt[j] > 30. && fabs(ohBJetL2CorrectedEta[j]) < 3.) { // ET and eta cuts
						
		    if (ohBJetIPL25Tag[j] >= 0)
		      { // Level 2.5 b tag  
			if (ohBJetIPL3Tag[j] >= 3.3)
			  { // Level 3 b tag  
			    rc++;
			  }
		      }
		  } // ET and eta cuts
		}//loop over jets
				
	      if (rc >= 1 && OpenHlt1MuonPassed(0., 14., 17., 2., 1, 2.1, 2.1, 1, 2)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /**********************************************/
	
  /*Electron-jet cross-triggers*/

  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdT_TrkIdT_DiJet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 2 && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation	 
											999.,	 
											999., // Track iso barrel, Track iso endcap	 
											999.,	 
											999., // Track/pT iso barrel, Track/pT iso endcap	 
											999.,	 
											999., // H/ET iso barrel, H/ET iso endcap	 
											999.,	 
											999., // E/ET iso barrel, E/ET iso endcap	 
											0.10,	 
											0.075, // H/E barrel, H/E endcap	 
											0.011,	 
											0.031, // cluster shape barrel, cluster shape endcap	 
											999.,	 
											999., // R9 barrel, R9 endcap	 
											0.008,	 
											0.008, // Deta barrel, Deta endcap	 
											0.07,	 
											0.05 // Dphi barrel, Dphi endcap	 
											)>=1)	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }	 
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdT_TrkIdT_TriJet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 3 && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation	 
											999.,	 
											999., // Track iso barrel, Track iso endcap	 
											999.,	 
											999., // Track/pT iso barrel, Track/pT iso endcap	 
											999.,	 
											999., // H/ET iso barrel, H/ET iso endcap	 
											999.,	 
											999., // E/ET iso barrel, E/ET iso endcap	 
											0.10,	 
											0.075, // H/E barrel, H/E endcap	 
											0.011,	 
											0.031, // cluster shape barrel, cluster shape endcap	 
											999.,	 
											999., // R9 barrel, R9 endcap	 
											0.008,	 
											0.008, // Deta barrel, Deta endcap	 
											0.07,	 
											0.05 // Dphi barrel, Dphi endcap	 
											)>=1)	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }	 
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdT_TrkIdT_QuadJet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 4 && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation	 
											999.,	 
											999., // Track iso barrel, Track iso endcap	 
											999.,	 
											999., // Track/pT iso barrel, Track/pT iso endcap	 
											999.,	 
											999., // H/ET iso barrel, H/ET iso endcap	 
											999.,	 
											999., // E/ET iso barrel, E/ET iso endcap	 
											0.10,	 
											0.075, // H/E barrel, H/E endcap	 
											0.011,	 
											0.031, // cluster shape barrel, cluster shape endcap	 
											999.,	 
											999., // R9 barrel, R9 endcap	 
											0.008,	 
											0.008, // Deta barrel, Deta endcap	 
											0.07,	 
											0.05 // Dphi barrel, Dphi endcap	 
											)>=1)	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }

  ///VBF Paths
  //------ 3E33 ------//
  
  //Ele27 - TighterEleIdIsol - CaloJetCor pT 35 25 - Deta 3 + Third Jet 20
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta3_Jet20") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
            if (OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedTriJetPassed(
					      35.,
					      25.,
                                              20.,
					      -1.,
					      true,
					      "Calo",
					      3.,
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }

  //Ele30 - TighterEleIdIsol - CaloJetCor pT 35 25 - Deta 3.5 + Third Jet 25 - backup
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele30_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta3p5_Jet25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
            if (OpenHlt1ElectronVbfEleIDPassed(30., 20., true, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedTriJetPassed(
					      35.,
					      25.,
                                              25.,
					      -1.,
					      true,
					      "Calo",
					      3.5,
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }

  //Ele27 - TighterEleId NO Isol - CaloJetCor pT 35 25 - Deta 3 + Third Jet 20 - background
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_TrkIdT_Jet35_Jet25_Deta3_Jet20") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
            if (OpenHlt1ElectronVbfEleIDPassed(27., 20., false, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedTriJetPassed(
					      35.,
					      25.,
                                              20.,
					      -1.,
					      true,
					      "Calo",
					      3.,
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }
  
  //------ 5E33 ------//
  
  //Ele27 - WP80 - PFJet pT 25 25 - Deta 3
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_WP80_DiPFJet25_Deta3") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, true)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
              25.,
              -1.,
              false,
              "PF",
              3.,
              false,
              false,
              EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }

  //Ele32 - WP80 - PFJet pT 25 25 - Deta 3p5
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele32_WP80_DiPFJet25_Deta3p5") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (OpenHlt1ElectronVbfEleIDPassed(32., 20., true, EtMaxIt, &EleIts, true)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
              25.,
              -1.,
              false,
              "PF",
              3.5,
              false,
              false,
              EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }
  
  //Ele27 - TighterEleId NO Isol - PFJet pT 25 25 - Deta 3 - background
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_TrkIdT_DiPFJet25_Deta3") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (OpenHlt1ElectronVbfEleIDPassed(27., 20., false, EtMaxIt, &EleIts, false)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
              25.,
              -1.,
              false,
              "PF",
              3.,
              false,
              false,
              EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }

  ///ggH Paths
  //------ 3E33 ------//

  //Ele27 - TighterEleIdIsol - CaloJetCor pT 30 25 - pfMHT 20
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
	    if (pfMHT > 20. &&
		OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedDiJetPassed(
					      30.,
					      25.,
					      2.6,
					      true,
					      "Calo",
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }
  
  //Ele30 - TighterEleIdIsol - CaloJetCor pT 30 30 - pfMHT 25 - backup
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele30_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30_PFMHT25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
	    if (pfMHT > 25. &&
		OpenHlt1ElectronVbfEleIDPassed(30., 20., true, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedDiJetPassed(
					      30.,
					      30.,
					      2.6,
					      true,
					      "Calo",
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }
   
  //Ele27 - TighterEleIdIsol - CaloJetCor pT 30 25 - pfMHT turn on studies
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
            if (OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedDiJetPassed(
					      30.,
					      25.,
					      2.6,
					      true,
					      "Calo",
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }

  //Ele27 - TighterEleId NO Isol - CaloJetCor pT 30 25 - background
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_TrkIdT_CentralJet30_CentralJet25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it) )
	  {
	    int EtMaxIt = -1;
	    std::vector<int> EleIts;
            if (OpenHlt1ElectronVbfEleIDPassed(27., 20., false, EtMaxIt, &EleIts, false)
		>= 1)
	      {
		if (OpenHltCleanedDiJetPassed(
					      30.,
					      25.,
					      2.6,
					      true,
					      "Calo",
					      0.,
					      false,
					      false,
					      EleIts) >= 1)
		  triggerBit[it] = true;
	      }
	  }
      }
  }
 
  
  //------ 5E33 ------//

  //Ele27 - WP80 - PFJet pT 25 25 - PFMHT15 - signal
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_WP80_DiCentralPFJet25_PFMHT15") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (pfMHT > 15. && 
            OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, true)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
          25.,
          2.6,
          false,
          "PF",
          0.,
          false,
          false,
          EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }
  
  //Ele32 - WP80 - PFJet pT 25 25 - PFMHT20 - backup
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele32_WP80_DiCentralPFJet25_PFMHT25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (pfMHT > 25. && 
            OpenHlt1ElectronVbfEleIDPassed(32., 20., true, EtMaxIt, &EleIts, true)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
          25.,
          2.6,
          false,
          "PF",
          0.,
          false,
          false,
          EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }

  
  //Ele27 - TighterEleId NO Isol - Central PFJet pT 25 25 - background
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_CaloIdVT_TrkIdT_DiCentralPFJet25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (OpenHlt1ElectronVbfEleIDPassed(27., 20., false, EtMaxIt, &EleIts, false)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              30.,
          25.,
          2.6,
          false,
          "PF",
          0.,
          false,
          false,
          EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }

  //Ele27 - WP80 - PFJet pT 25 25 - PFMHT turn on
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_Ele27_WP80_DiCentralPFJet25") == 0) {
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second == 1)
    {
      if (prescaleResponse(menu, cfg, rcounter, it) )
      {
        int EtMaxIt = -1;
        std::vector<int> EleIts;
        if (OpenHlt1ElectronVbfEleIDPassed(27., 20., true, EtMaxIt, &EleIts, true)
            >= 1)
        {
          if (OpenHltCleanedDiJetPassed(
              25.,
          25.,
          2.6,
          false,
          "PF",
          0.,
          false,
          false,
          EleIts) >= 1)
            triggerBit[it] = true;
        }
      }
    }
  }
  
  // 2011-05-11: promoted to v4, removed R9 cut, TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation  
						  999.,
						  999., // Track iso barrel, Track iso endcap  
						  999,
						  999, // Track/pT iso barrel, Track/pT iso endcap  
						  0.2,
						  0.2, // H/ET iso barrel, H/ET iso endcap  
						  0.2,
						  0.2, // E/ET iso barrel, E/ET iso endcap  
						  0.15,
						  0.1, // H/E barrel, H/E endcap  
						  0.014,
						  0.035, // cluster shape barrel, cluster shape endcap  
						  999.,//0.98,
						  999.,//1.0, // R9 barrel, R9 endcap  
						  999,
						  999, // Deta barrel, Deta endcap  
						  999,
						  999 // Dphi barrel, Dphi endcap  
						  )>=1 && OpenHlt1CorJetPassed(40)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

/*
HLT_Photon70_CaloIdL_HT300:
	Pt 70 naturally
	shape filter: 0.014, 0.035
	H/E filter: 0.15, 0.1
	min Ht 300; minPtJet 40; etaJet 3.0; mode 4
HLT_Photon70_CaloIdL_MHT70:
        Pt 70 naturally
        shape filter: 0.014, 0.035
        H/E filter: 0.15, 0.1
        minMht 300; minPtJet 30; etaJet 3.0; mode 1	
*/

//2011-06-08: These should be updated once pfMHTphi becmoes available in ntuple
else if (triggerName.CompareTo("Ele25_CaloIdWP80_CaloIsoWP80_TrkIdWP80_TrkIsoWP80_PFMT40") == 0 || 
	triggerName.CompareTo("OpenHLT_Ele25_WP80_PFMT40") == 0 )//two names, one definition for historical reasons. 
{
	if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
		if (prescaleResponse(menu, cfg, rcounter, it))
		{
			std::vector<int> eleVector = VectorOpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  x
											   999.,
											   999., // Track iso barrel, Track iso endcap  
											   0.09,
											   0.04, // Track/pT iso barrel, Track/pT iso endcap  x
											   0.1,
											   0.025, // H/ET iso barrel, H/ET iso endcap  x
											   0.07,
											   0.05, // E/ET iso barrel, E/ET iso endcap x
											   0.04,
											   0.025, // H/E barrel, H/E endcap  x
											   0.01,
											   0.03, // cluster shape barrel, cluster shape endcap x
											   999.,
											   999., // R9 barrel, R9 endcap  
											   0.004,
											   0.007, // Deta barrel, Deta endcap  x
											   0.06,
											   0.03 // Dphi barrel, Dphi endcap  x
											   );
			if(eleVector.size()>=1)
			{
					//calculate pfMHTphi. Once the new producer module is in place, section should be replaced by the pfMHTPhi in ntuple
				float pfMHTphi = 0.;
				TLorentzVector sum;
				TLorentzVector thisjet;
				for(int ijet = 0; ijet<NohPFJet; ijet++){
					thisjet.SetPtEtaPhiM(pfJetPt[ijet], pfJetEta[ijet], pfJetPhi[ijet], 0.);
					sum += thisjet;
				}
				pfMHTphi = (-sum).Phi();
					//end calculate pfMHTphi

					
				//do any of the passing ele make pfMT 40 with pfMHT?
				for (unsigned int iele = 0; iele<eleVector.size(); iele++) {
					int i = eleVector[iele];
					float dphi = fabs(pfMHTphi-ohElePhi[i]);
					if (dphi > 3.14159){ dphi = (2.0 * 3.14159) - dphi; }
					float pfMT = sqrt(2.*ohEleEt[i]*pfMHT*(1.-cos(dphi)));
					if (pfMT > 40) {
						triggerBit[it] = true;
						break;
					}
				}//end for all ele
			}//end if there are passing electrons
		}
	}
}

//2011-06-08: These should be updated once pfMHTphi becmoes available in ntuple
else if (triggerName.CompareTo("Ele27_CaloIdWP70_CaloIsoWP70_TrkIdWP70_TrkIsoWP70_PFMT40_PFMHT20")== 0 || 
	triggerName.CompareTo("OpenHLT_Ele27_WP70_PFMT40_PFMHT20")== 0 ) //two names, one definition for historical reasons. 
{
	if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
		if (prescaleResponse(menu, cfg, rcounter, it))
		{
			
			if (pfMHT > 20.){
				std::vector<int> eleVector = VectorOpenHlt1ElectronSamHarperPassed(
												   27., 0, // ET, L1isolation 
												   999.,
												   999., // Track iso barrel, Track iso endcap  
												   0.05,
												   0.025, // Track/pT iso barrel, Track/pT iso endcap
												   0.03,
												   0.02, // H/ET iso barrel, H/ET iso endcap
												   0.06,
												   0.025, // E/ET iso barrel, E/ET iso endcap
												   0.025,
												   0.025, // H/E barrel, H/E endcap
												   0.01,
												   0.03, // cluster shape barrel, cluster shape endcap  x
												   999.,
												   999., // R9 barrel, R9 endcap  
												   0.004,
												   0.005, // Deta barrel, Deta endcap
												   0.03,
												   0.02 // Dphi barrel, Dphi endcap
												   );
				if(eleVector.size()>=1)
				{
						//calculate pfMHTphi. Once the new producer module is in place, section should be replaced by the pfMHTPhi in ntuple
					float pfMHTphi = 0.;
					TLorentzVector sum;
					TLorentzVector thisjet;
					for(int ijet = 0; ijet<NohPFJet; ijet++){
						thisjet.SetPtEtaPhiM(pfJetPt[ijet], pfJetEta[ijet], pfJetPhi[ijet], 0.);
						sum += thisjet;
					}
					pfMHTphi = (-sum).Phi();
						//end calculate pfMHTphi
					
					
					//do any of the passing ele make pfMT 40 with pfMHT?
					for (unsigned int iele = 0; iele<eleVector.size(); iele++) { 
						int i = eleVector[iele];
						float dphi = fabs(pfMHTphi-ohElePhi[i]);
						if (dphi > 3.14159){ dphi = (2.0 * 3.14159) - dphi; }
						float pfMT = sqrt(2.*ohEleEt[i]*pfMHT*(1.-cos(dphi))); 
						if (pfMT > 40) {
							triggerBit[it] = true;
							break;
						}
					}//end for all ele
				}//end if there are passing electrons
			}//end pfMHT cut
		}
	}
}

//2011-06-08: These should be updated once pfMHTphi becmoes available in ntuple
else if (triggerName.CompareTo("OpenHLT_Ele27_WP80_PFMT50_v1") == 0)
{
	if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
		if (prescaleResponse(menu, cfg, rcounter, it))
		{
			std::vector<int> eleVector = VectorOpenHlt1ElectronSamHarperPassed(27., 0, // ET, L1isolation  x
											   999.,
											   999., // Track iso barrel, Track iso endcap  
											   0.09,
											   0.04, // Track/pT iso barrel, Track/pT iso endcap  x
											   0.1,
											   0.025, // H/ET iso barrel, H/ET iso endcap  x
											   0.07,
											   0.05, // E/ET iso barrel, E/ET iso endcap x
											   0.04,
											   0.025, // H/E barrel, H/E endcap  x
											   0.01,
											   0.03, // cluster shape barrel, cluster shape endcap x
											   999.,
											   999., // R9 barrel, R9 endcap  
											   0.004,
											   0.007, // Deta barrel, Deta endcap  x
											   0.06,
											   0.03 // Dphi barrel, Dphi endcap  x
											   );
			if(eleVector.size()>=1)
			{
					//calculate pfMHTphi. Once the new producer module is in place, section should be replaced by the pfMHTPhi in ntuple
				float pfMHTphi = 0.;
				TLorentzVector sum;
				TLorentzVector thisjet;
				for(int ijet = 0; ijet<NohPFJet; ijet++){
					thisjet.SetPtEtaPhiM(pfJetPt[ijet], pfJetEta[ijet], pfJetPhi[ijet], 0.);
					sum += thisjet;
				}
				pfMHTphi = (-sum).Phi();
					//end calculate pfMHTphi
				
					//do any of the passing ele make pfMT 40 with pfMHT?
				for (unsigned int iele = 0; iele<eleVector.size(); iele++) {
					int i = eleVector[iele];
					float dphi = fabs(pfMHTphi-ohElePhi[i]);
					if (dphi > 3.14159){ dphi = (2.0 * 3.14159) - dphi; }
					float pfMT = sqrt(2.*ohEleEt[i]*pfMHT*(1.-cos(dphi)));
					if (pfMT > 50) {
						triggerBit[it] = true;
						break;
					}
				}//end for all ele
			}//end if there are passing electrons
		}
	}
}

//2011-06-08: These should be updated once pfMHTphi becmoes available in ntuple
else if (triggerName.CompareTo("OpenHLT_Ele32_WP70_PFMT50_v1")  == 0)
{
	if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
		if (prescaleResponse(menu, cfg, rcounter, it)) 
		{
			
			if (pfMHT > 20.){
				std::vector<int> eleVector = VectorOpenHlt1ElectronSamHarperPassed(
												   32., 0, // ET, L1isolation 
												   999.,
												   999., // Track iso barrel, Track iso endcap  
												   0.05,
												   0.025, // Track/pT iso barrel, Track/pT iso endcap
												   0.03,
												   0.02, // H/ET iso barrel, H/ET iso endcap
												   0.06,
												   0.025, // E/ET iso barrel, E/ET iso endcap
												   0.025,
												   0.025, // H/E barrel, H/E endcap
												   0.01,
												   0.03, // cluster shape barrel, cluster shape endcap  x
												   999.,
												   999., // R9 barrel, R9 endcap  
												   0.004,
												   0.005, // Deta barrel, Deta endcap
												   0.03,
												   0.02 // Dphi barrel, Dphi endcap
												   );
				if(eleVector.size()>=1)
				{
						//calculate pfMHTphi. Once the new producer module is in place, section should be replaced by the pfMHTPhi in ntuple
					float pfMHTphi = 0.;
					TLorentzVector sum;
					TLorentzVector thisjet;
					for(int ijet = 0; ijet<NohPFJet; ijet++){
						thisjet.SetPtEtaPhiM(pfJetPt[ijet], pfJetEta[ijet], pfJetPhi[ijet], 0.);
						sum += thisjet;
					}
					pfMHTphi = (-sum).Phi();
						//end calculate pfMHTphi
					
						//do any of the passing ele make pfMT 40 with pfMHT?
					for (unsigned int iele = 0; iele<eleVector.size(); iele++) { 
						int i = eleVector[iele];
						float dphi = fabs(pfMHTphi-ohElePhi[i]);
						if (dphi > 3.14159){ dphi = (2.0 * 3.14159) - dphi; }
						float pfMT = sqrt(2.*ohEleEt[i]*pfMHT*(1.-cos(dphi)));
						if (pfMT > 50) {
							triggerBit[it] = true;
							break;
						}
					}//end for all ele
				}//end if there are passing electrons
			}//end pfMHT cut
		}
	}
}

	
  /**********************************************/
	
  // Single Top Triggers : Electron channel
	
  /**********************************************/
  else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation 
						  999.,
						  999., // Track iso barrel, Track iso endcap 
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap 
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap 
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap 
						  0.05,
						  0.05, // H/E barrel, H/E endcap 
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap 
						  0.98,
						  1.0, // R9 barrel, R9 endcap 
						  0.008,
						  0.008, // Deta barrel, Deta endcap 
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap 
						  )>=1 && OpenHlt1BJetPassedEleRemoval(30., 3.0, 0.3, // jet ET, eta, DrCut
										       0.,
										       3.3, // discL25, discL3
										       25.,
										       0, // ET, L1isolation 
										       999.,
										       999., // Track iso barrel, Track iso endcap 
										       999.,
										       999., // Track/pT iso barrel, Track/pT iso endcap 
										       999.,
										       999., // H/ET iso barrel, H/ET iso endcap 
										       999.,
										       999., // E/ET iso barrel, E/ET iso endcap 
										       0.05,
										       0.05, // H/E barrel, H/E endcap 
										       0.011,
										       0.031, // cluster shape barrel, cluster shape endcap 
										       0.98,
										       1.0, // R9 barrel, R9 endcap 
										       0.008,
										       0.008, // Deta barrel, Deta endcap 
										       0.07,
										       0.05 // Dphi barrel, Dphi endcap 
										       )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  //Isolated Ele path
  else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation 
						  999.,
						  999., // Track iso barrel, Track iso endcap 
						  0.125,
						  0.075, // Track/pT iso barrel, Track/pT iso endcap 
						  0.125,
						  0.075, // H/ET iso barrel, H/ET iso endcap 
						  0.125,
						  0.075, // E/ET iso barrel, E/ET iso endcap 
						  0.05,
						  0.05, // H/E barrel, H/E endcap 
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap 
						  0.98,
						  1.0, // R9 barrel, R9 endcap 
						  0.008,
						  0.008, // Deta barrel, Deta endcap 
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap 
						  )>=1 && OpenHlt1BJetPassedEleRemoval(30., 3.0, 0.3, // jet ET, eta, DrCut
										       0.,
										       3.3, // discL25, discL3
										       25.,
										       0, // ET, L1isolation 
										       999.,
										       999., // Track iso barrel, Track iso endcap 
										       0.125,
										       0.075, // Track/pT iso barrel, Track/pT iso endcap 
										       0.125,
										       0.075, // H/ET iso barrel, H/ET iso endcap 
										       0.125,
										       0.075, // E/ET iso barrel, E/ET iso endcap 
										       0.05,
										       0.05, // H/E barrel, H/E endcap 
										       0.011,
										       0.031, // cluster shape barrel, cluster shape endcap 
										       0.98,
										       1.0, // R9 barrel, R9 endcap 
										       0.008,
										       0.008, // Deta barrel, Deta endcap 
										       0.07,
										       0.05 // Dphi barrel, Dphi endcap
										       )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation 
						  999.,
						  999., // Track iso barrel, Track iso endcap 
						  999.,
						  999., // Track/pT iso barrel, Track/pT iso endcap 
						  999.,
						  999., // H/ET iso barrel, H/ET iso endcap 
						  999.,
						  999., // E/ET iso barrel, E/ET iso endcap 
						  0.05,
						  0.05, // H/E barrel, H/E endcap 
						  0.011,
						  0.031, // cluster shape barrel, cluster shape endcap 
						  0.98,
						  1.0, // R9 barrel, R9 endcap 
						  0.008,
						  0.008, // Deta barrel, Deta endcap 
						  0.07,
						  0.05 // Dphi barrel, Dphi endcap 
						  )>=1 && OpenHlt1BJetPassedEleRemoval(30., 3.0, 0.3, // jet ET, eta, DrCut
										       0.,
										       2.0, // discL25, discL3
										       25.,
										       0, // ET, L1isolation 
										       999.,
										       999., // Track iso barrel, Track iso endcap 
										       999.,
										       999., // Track/pT iso barrel, Track/pT iso endcap 
										       999.,
										       999., // H/ET iso barrel, H/ET iso endcap 
										       999.,
										       999., // E/ET iso barrel, E/ET iso endcap 
										       0.05,
										       0.05, // H/E barrel, H/E endcap 
										       0.011,
										       0.031, // cluster shape barrel, cluster shape endcap 
										       0.98,
										       1.0, // R9 barrel, R9 endcap 
										       0.008,
										       0.008, // Deta barrel, Deta endcap 
										       0.07,
										       0.05 // Dphi barrel, Dphi endcap 
										       )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

   /* e + 1/2/3/4 jets */

   // 2011-05-03 fixed the jets cleaning according to ConfDB implementation
   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_CentralJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNCorJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
  	          0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_DiCentralJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNCorJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNCorJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_QuadCentralJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNCorJetPassedEleRemoval(4, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   // ele_PFjets triggers

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_CentralPFJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNPFJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
  	          0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_DiCentralPFJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNPFJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_TriCentralPFJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNPFJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_TrkIdT_QuadCentralPFJet30") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation  
                  999.,
                  999., // Track iso barrel, Track iso endcap  
                  999.,
                  999., // Track/pT iso barrel, Track/pT iso endcap  
                  999.,
                  999., // H/ET iso barrel, H/ET iso endcap  
                  999.,
                  999., // E/ET iso barrel, E/ET iso endcap  
                  0.05,
                  0.05, // H/E barrel, H/E endcap  
                  0.011,
                  0.031, // cluster shape barrel, cluster shape endcap  
                  0.98,
                  1.0, // R9 barrel, R9 endcap  
                  0.008,
                  0.008, // Deta barrel, Deta endcap  
                  0.07,
                  0.05 // Dphi barrel, Dphi endcap  
             )>=1 && OpenHltNPFJetPassedEleRemoval(4, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
  		  25.,
                  0, // ET, L1isolation
   		  999.,
		  999., // Track iso barrel, Track iso endcap  
		  999.,
		  999., // Track/pT iso barrel, Track/pT iso endcap  
		  999.,
		  999., // H/ET iso barrel, H/ET iso endcap  
		  999.,
		  999., // E/ET iso barrel, E/ET iso endcap  
		  0.05,
		  0.05, // H/E barrel, H/E endcap  
		  0.011,
		  0.031, // cluster shape barrel, cluster shape endcap  
		  0.98,
		  1.0, // R9 barrel, R9 endcap  
		  0.008,
		  0.008, // Deta barrel, Deta endcap  
		  0.07,
		  0.05 // Dphi barrel, Dphi endcap
	    ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

  ///////////

  //May5th
 else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
						)>=1 )
            {
               triggerBit[it] = true;
            }
         }
      }
   }



   /* loose iso/ loose ID e + 1/2/3 jets */
   //added for 1E33 to be deployed on 2011-05-19

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_CentralJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_CentralDiJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_CentralTriJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_CentralQuadJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(4, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.2,
                  0.2, // Track/pT iso barrel, Track/pT iso endcap
                  0.2,
                  0.2, // H/ET iso barrel, H/ET iso endcap
                  0.2,
                  0.2, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }

   /* tight iso/ loose ID e + 1/2/3/4 jets */
   //added for 1E33 to be deployed on 2011-05-19

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoT_TrkIdVL_TrkIsoT_CentralJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }
      
   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoT_TrkIdVL_TrkIsoT_CentralDiJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }
   
   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoT_TrkIdVL_TrkIsoT_CentralTriJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }
   
   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdL_CaloIsoT_TrkIdVL_TrkIsoT_CentralQuadJet30_v3") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
             )>=1 && OpenHltNCorJetPassedEleRemoval(4, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
                  25.,
                  0, // ET, L1isolation
                  999.,
                  999., // Track iso barrel, Track iso endcap
                  0.125,
                  0.075, // Track/pT iso barrel, Track/pT iso endcap
                  0.125,
                  0.075, // H/ET iso barrel, H/ET iso endcap
                  0.125,
                  0.075, // E/ET iso barrel, E/ET iso endcap
                  0.15,
                  0.10, // H/E barrel, H/E endcap
                  0.014,
                  0.035, // cluster shape barrel, cluster shape endcap
                  0.98,
                  1.0, // R9 barrel, R9 endcap
                  0.01,
                  0.01, // Deta barrel, Deta endcap
                  0.15,
                  0.10 // Dphi barrel, Dphi endcap
            ))
            {
               triggerBit[it] = true;
            }
         }
      }
   }
   
  /* tight iso/ tight ID e + 1/2/3/4 jets */
  //added for 1E33 to be deployed on 2011-05-19

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNCorJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNCorJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNCorJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_QuadCentralJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNCorJetPassedEleRemoval(4, 30., 2.6, 0.3, //jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }


  //PF Jets
  /* tight iso/ tight ID e + 1/2/3/4 jets */

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNPFJetPassedEleRemoval(1, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralPFJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNPFJetPassedEleRemoval(2, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNPFJetPassedEleRemoval(3, 30., 2.6, 0.3, //N, jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }

   else if (triggerName.CompareTo("OpenHLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_QuadCentralPFJet30") == 0)
     {
       if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	 {
	   if (prescaleResponse(menu, cfg, rcounter, it))
	     {
	       if (OpenHlt1ElectronSamHarperPassed(25., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   0.125,
						   0.075, // Track/pT iso barrel, Track/pT iso endcap
						   0.125,
						   0.075, // H/ET iso barrel, H/ET iso endcap
						   0.125,
						   0.075, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap
						   999.,
						   999., // R9 barrel, R9 endcap
						   0.008,
						   0.008, // Deta barrel, Deta endcap
						   0.07,
						   0.05 // Dphi barrel, Dphi endcap
						   )>=1 && OpenHltNPFJetPassedEleRemoval(4, 30., 2.6, 0.3, //jetPt, jetEta, dR cut
											  25.,
											  0, // ET, L1isolation
											  999.,
											  999., // Track iso barrel, Track iso endcap
											  0.125,
											  0.075, // Track/pT iso barrel, Track/pT iso endcap
											  0.125,
											  0.075, // H/ET iso barrel, H/ET iso endcap
											  0.125,
											  0.075, // E/ET iso barrel, E/ET iso endcap
											  0.05,
											  0.05, // H/E barrel, H/E endcap
											  0.011,
											  0.031, // cluster shape barrel, cluster shape endcap
											  999.,
											  999., // R9 barrel, R9 endcap
											  0.008,
											  0.008, // Deta barrel, Deta endcap
											  0.07,
											  0.05 // Dphi barrel, Dphi endcap
											  ))
		 {
		   triggerBit[it] = true;
		 }
	     }
	 }
     }



	
  /* Minbias */
  else if (triggerName.CompareTo("OpenHLT_MinBiasBSC_OR") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      bool techTriggerBSCOR = (bool) L1Tech_BSC_minBias_OR_v0;
	      if (techTriggerBSCOR)
		triggerBit[it] = true;
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_MinBiasBSC") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      bool techTriggerBSC1 = (bool) L1Tech_BSC_minBias_threshold1_v0;
	      bool techTriggerBSC2 = (bool) L1Tech_BSC_minBias_threshold2_v0;
	      bool techTriggerBS3 = (bool) L1Tech_BSC_minBias_inner_threshold1_v0;
	      bool techTriggerBS4 = (bool) L1Tech_BSC_minBias_inner_threshold2_v0;
				
	      if (techTriggerBSC1 || techTriggerBSC2 || techTriggerBS3
		  || techTriggerBS4)
		triggerBit[it] = true;
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_MinBiasPixel_SingleTrack") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PixelTrackPassed(0.0, 1.0, 0.0)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_MinBiasPixel_DoubleTrack") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PixelTrackPassed(0.0, 1.0, 0.0)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_MinBiasPixel_DoubleIsoTrack5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1PixelTrackPassed(5.0, 1.0, 1.0)>=0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_ZeroBias") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (map_BitOfStandardHLTPath.find("OpenL1_ZeroBias")->second == 1)
		triggerBit[it] = true;
	    }
	}
    }
	
  /* AlCa */
  else if (triggerName.CompareTo("OpenAlCa_HcalPhiSym") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (ohHighestEnergyHFRecHit > 0 || ohHighestEnergyHBHERecHit > 0)
		{
		  // Require one RecHit with E > 0 MeV in any HCAL subdetector
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoTrackHB") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      bool passL2=false;
	      for (int itrk=0; itrk<ohIsoPixelTrackHBL2N; itrk++)
		{
		  if (ohIsoPixelTrackHBL2P[itrk]>8.0
		      && TMath::Abs(ohIsoPixelTrackHBL2Eta[itrk])>0.0
		      && TMath::Abs(ohIsoPixelTrackHBL2Eta[itrk])<1.3
		      && ohIsoPixelTrackHBL2MaxNearP[itrk]<2.0)
		    passL2=true;
		}
				
	      bool passL3=false;
	      for (int itrk=0; itrk<ohIsoPixelTrackHBL3N; itrk++)
		{
		  if (ohIsoPixelTrackHBL3P[itrk]>20.0
		      && TMath::Abs(ohIsoPixelTrackHBL3Eta[itrk])>0.0
		      && TMath::Abs(ohIsoPixelTrackHBL3Eta[itrk])<1.3
		      && ohIsoPixelTrackHBL3MaxNearP[itrk]<2.0)
		    passL3=true;
		}
				
	      if (passL2 && passL3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoTrackHE") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
				
	      bool passL2=false;
	      for (int itrk=0; itrk<ohIsoPixelTrackHEL2N; itrk++)
		{
		  if (ohIsoPixelTrackHEL2P[itrk]>12.0
		      && TMath::Abs(ohIsoPixelTrackHEL2Eta[itrk])>0.0
		      && TMath::Abs(ohIsoPixelTrackHEL2Eta[itrk])<2.2
		      && ohIsoPixelTrackHEL2MaxNearP[itrk]<2.0)
		    passL2=true;
		}
				
	      bool passL3=false;
	      for (int itrk=0; itrk<ohIsoPixelTrackHEL3N; itrk++)
		{
		  if (ohIsoPixelTrackHEL3P[itrk]>38.0
		      && TMath::Abs(ohIsoPixelTrackHEL3Eta[itrk])>0.0
		      && TMath::Abs(ohIsoPixelTrackHEL3Eta[itrk])<2.2
		      && ohIsoPixelTrackHEL3MaxNearP[itrk]<2.0)
		    passL3=true;
		}
				
	      if (passL2 && passL3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*Photon-electron cross-triggers*/
  // Not finished yet. 
	
	
  // 2011-03-29: promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1ElectronSamHarperPassed(
										   8., 0, // ET, L1isolation  
										   999.,
										   999., // Track iso barrel, Track iso endcap  
										   999,
										   999, // Track/pT iso barrel, Track/pT iso endcap  
										   0.2,
										   0.2, // H/ET iso barrel, H/ET iso endcap  
										   0.2,
										   0.2, // E/ET iso barrel, E/ET iso endcap  
										   0.15,
										   0.1, // H/E barrel, H/E endcap  
										   0.014,
										   0.035, // cluster shape barrel, cluster shape endcap  
										   0.98,
										   1.0, // R9 barrel, R9 endcap  
										   999,
										   999, // Deta barrel, Deta endcap  
										   999,
										   999); // Dphi barrel, Dphi endcap  
	      if (firstVector.size()>=1){
		std::vector<int> secondVector = VectorOpenHlt1PhotonSamHarperPassed(
										    20., 0, // ET, L1isolation
										    999.,
										    999., // Track iso barrel, Track iso endcap
										    3.0,
										    3.0, // Track/pT iso barrel, Track/pT iso endcap
										    3.0,
										    3.0, // H/ET iso barrel, H/ET iso endcap
										    5.0,
										    5.0, // E/ET iso barrel, E/ET iso endcap
										    0.05,
										    0.05, // H/E barrel, H/E endcap  
										    0.011,
										    0.031, // cluster shape barrel, cluster shape endcap  
										    0.98,
										    999., // R9 barrel, R9 endcap
										    999.,
										    999., // Deta barrel, Deta endcap
										    999.,
										    999.); // Dphi barrel, Dphi endcap
		if (secondVector.size()>=1){
                        
		  // make sure they're not the same object
		  TLorentzVector ele;
		  TLorentzVector pho;
		  float deltaR;
		  int NMatch = 0;
		  for (unsigned int i=0; i<firstVector.size(); i++)
		    {
		      ele.SetPtEtaPhiM(
				       ohEleEt[firstVector[i]],
				       ohEleEta[firstVector[i]],
				       ohElePhi[firstVector[i]],
				       0.);
		      for (unsigned int j=0; j<secondVector.size() ; j++)
			{
			  pho.SetPtEtaPhiM(
					   ohPhotEt[secondVector[j]],
					   ohPhotEta[secondVector[j]],
					   ohPhotPhi[secondVector[j]],
					   0.);
			  deltaR = ele.DeltaR(pho);
			  if (deltaR < 0.3)
			    {
			      NMatch++;
			    }
			}
		    }			    
						
		  if (((firstVector.size() -  NMatch) != 0) && ((secondVector.size() - NMatch ) != 0))
		    {
		      triggerBit[it] = true;
		    }
						
		}
	      }
	    }
	}
    }
	
  /* muon-jet/MET/HT cross-triggers */
  else if (triggerName.CompareTo("OpenHLT_Mu3_Dijet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 2 && OpenHlt1MuonPassed(0., 0., 3., 2., 0)>=1 )	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }	 
  else if (triggerName.CompareTo("OpenHLT_Mu3_Trijet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 3 && OpenHlt1MuonPassed(0., 0., 3., 2., 0)>=1 )	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }	 
  else if (triggerName.CompareTo("OpenHLT_Mu3_Quadjet30") == 0)	 
    {	 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)	 
	{	 
	  if (prescaleResponse(menu, cfg, rcounter, it))	 
	    {	 
	      if (OpenHlt1CorJetPassed(30., 3.) >= 4 && OpenHlt1MuonPassed(0., 0., 3., 2., 0)>=1 )	 
		{	 
		  triggerBit[it] = true;	 
		}	 
	    }	 
	}	 
    }	 

  else if (triggerName.CompareTo("OpenHLT_Mu17_CentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 17., 2., 0)>=1
		  && OpenHlt1CorJetPassed( 30, 2.6)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu17_DiCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 17., 2., 0)>=1
		  && OpenHlt1CorJetPassed( 30, 2.6)>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu17_TriCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 17., 2., 0)>=1
		  && OpenHlt1CorJetPassed( 30, 2.6)>=3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu17_QuadCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(10., 10., 17., 2., 0)>=1
                  && OpenHlt1CorJetPassed( 30, 2.6)>=4)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_CentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(0., 10., 17., 2., 1)>=1
                  && OpenHlt1CorJetPassed( 30, 2.6)>=1)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_DiCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(0., 10., 17., 2., 1)>=1
                  && OpenHlt1CorJetPassed( 30, 2.6)>=2)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu20_DiCentralJet34") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(0., 10., 20., 2., 1)>=1
                  && OpenHlt1CorJetPassed( 34, 2.6)>=2)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_TriCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(0., 10., 17., 2., 1)>=1
                  && OpenHlt1CorJetPassed( 30, 2.6)>=3)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_QuadCentralJet30") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(0., 10., 17., 2., 1)>=1
                  && OpenHlt1CorJetPassed( 30, 2.6)>=4)
                {
                  triggerBit[it] = true;
                }
            }
        }
    }


// ttbar lep+jets PFjet paths
  else if (triggerName.CompareTo("OpenHLT_Mu17_eta2p1_CentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 0, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=1)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_Mu17_eta2p1_DiCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 0, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=2)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_Mu17_eta2p1_TriCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 0, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=3)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_Mu17_eta2p1_QuadCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 0, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=4)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 

  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1_CentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 1, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=1)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1_DiCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 1, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=2)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1_TriCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 1, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=3)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 
  else if (triggerName.CompareTo("OpenHLT_IsoMu17_eta2p1_QuadCentralPFJet30") == 0) 
    { 
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1) 
        { 
          if (prescaleResponse(menu, cfg, rcounter, it)) 
            { 
	      if (OpenHlt1MuonPassed(0., 14., 17., 2., 1, 2.1, 2.1, 1, 2)>=1 
                  && OpenHlt1PFJetPassed(30, 2.6)>=4)
                { 
                  triggerBit[it] = true; 
                } 
            } 
        } 
    } 

  else if (isDoubleMuX_HTXTrigger(triggerName, thresholds))
     {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., thresholds[0], 2., 0)>=2
		  && OpenHltSumCorHTPassed( thresholds[1])>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDoubleMuX_MassX_HTXTrigger(triggerName, thresholds))
     {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      int nMu = OpenHlt1MuonPassed(0., 0., thresholds[0], 2., 0);
	      if (nMu >=2
		  && OpenHltInvMassCutMu(nMu, thresholds[1])
		  && OpenHltSumCorHTPassed( thresholds[2])>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

 else if (isDoubleMuX_MassX_HTFJX_PFHTXTrigger(triggerName, thresholds))
     {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      int nMu = OpenHlt1MuonPassed(0., 0., thresholds[0], 2., 0);
	      if (nMu >=2
		  && OpenHltInvMassCutMu(nMu, thresholds[1])
		  && OpenHltSumFJCorHTPassed( thresholds[2])>0
		  && OpenHltSumPFHTPassed( thresholds[3])>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDoubleEleX_MassX_HTXTrigger(triggerName,  caloId,  caloIso,  trkId,  trkIso, thresholds))
     {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      int nEle = OpenHlt1ElectronPassed(thresholds[0], 
						 map_EGammaCaloId[caloId[0]],
						 map_EleCaloIso[caloIso[0]],
						 map_EleTrkId[trkId[0]],
						 map_EleTrkIso[trkIso[0]]
						 );
	      if (nEle >=2
		  && OpenHltInvMassCutEle(nEle, thresholds[1])
		  && OpenHltSumCorHTPassed( thresholds[2])>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	

  else if (isMuX_EleX_MassX_HTXTrigger(triggerName, caloId,  caloIso,  trkId,  trkIso, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      // int nMu = OpenHlt1MuonPassed(map_muThresholds[thresholds[0]],2.,0);
	      int nMu = OpenHlt1MuonPassed(0., 0., thresholds[0],2.,0);
	      int nEle = OpenHlt1ElectronPassed(thresholds[1], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
						);
	    
	      if (nMu >= 1 &&  nEle >= 1 && OpenHltInvMassCutEleMu(nEle, nMu, thresholds[2])
		  && OpenHltSumCorHTPassed(thresholds[3])>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
     }


  //AGB - HT + single electron + MET
	
  else if (isHTX_EleX_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_pfMHTXTrigger(triggerName, thresholds)) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(thresholds[2]) && OpenHltSumCorHTPassed(thresholds[0],40.,3.)>0  && (OpenHlt1ElectronSamHarperPassed(thresholds[1],0,          // ET, L1isolation
															     999., 999.,       // Track iso barrel, Track iso endcap
															     0.2, 0.2,        // Track/pT iso barrel, Track/pT iso endcap
															     0.2, 0.2,       // H/ET iso barrel, H/ET iso endcap
															     0.2, 0.2,       // E/ET iso barrel, E/ET iso endcap
															     0.15, 0.10,       // H/E barrel, H/E endcap
															     0.024, 0.040,       // cluster shape barrel, cluster shape endcap
															     0.98, 1.0,       // R9 barrel, R9 endcap
															     0.01, 0.01,       // Deta barrel, Deta endcap
															     0.15, 0.10        // Dphi barrel, Dphi endcap
															     )>=1))
	  {
	    triggerBit[it] = true;
	  }  
      }
    }
  }
	
	
	
	
  //AGB - HT + mu (+MET)
	
  else if (isHTX_MuX_pfMHTXTrigger(triggerName, thresholds)) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(thresholds[2]) && OpenHltSumCorHTPassed(thresholds[0])>0 && OpenHlt1MuonPassed(map_muThresholds[thresholds[1]],2.,0)>0)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }
	
 
  else if (isMuX_HTXTrigger(triggerName, thresholds)) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if( OpenHltSumCorHTPassed(thresholds[1])>0 && OpenHlt1MuonPassed(map_muThresholds[thresholds[0]],2.,0)>0)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }
	
    else if (isNJetPtTrigger(triggerName, thresholds)) {
        if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1) {
            if (prescaleResponse(menu,cfg,rcounter,it)) {
                int N= int(thresholds[0]);
                double pt= thresholds[1];
                if(OpenHltNJetPtPassed(N, pt)) {
                    //cout << triggerName << " (" << N << ", " << pt << ") passed" << endl;
                    triggerBit[it] = true;
                }
            }
        }
    }

    else if (isNTowCalEt0pTrigger(triggerName, thresholds)) {
        if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1) {
            if (prescaleResponse(menu,cfg,rcounter,it)) {
                int N= int(thresholds[0]);
                double Et= thresholds[1];
                if(OpenHltNTowCalEtPassed(N, Et)) {
                    //cout << triggerName << " (" << N << ", " << Et << ") passed" << endl;
                    triggerBit[it] = true;
                }
            }
        }
    }
 
  //AGB - HT + 2Tau (+MET)
	
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_HT200_DoubleLooseIsoPFTau10trk3_PFMHT35") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(35.) && 
	   OpenHltSumCorHTPassed(200.,40.,3.)>0 && 
	   OpenHltPFTauPassedNoMuonIDNoEleID(10.,3.,1.5,999.,
					     3.,3.,5.,2.,0, //muonID
					     5.,0,          // ET, L1isolation
					     999., 999.,       // Track iso barrel, Track iso endcap
					     999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
					     999., 999.,       // H/ET iso barrel, H/ET iso endcap
					     999., 999.,       // E/ET iso barrel, E/ET iso endcap
					     999., 999.,       // H/E barrel, H/E endcap
					     999., 999.,       // cluster shape barrel, cluster shape endcap
					     0.98, 1.0,       // R9 barrel, R9 endcap
					     999., 999.,       // Deta barrel, Deta endcap
					     999., 999.        // Dphi barrel, Dphi endcap
					     )>1)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }

  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_HT200_DoubleTightIsoPFTau10trk3_PFMHT35") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(35.) && 
	   OpenHltSumCorHTPassed(200.,40.,3.)>0 && 
	   OpenHltPFTauPassedNoMuonIDNoEleID(10.,3.,1.,1.5,
					     3.,3.,5.,2.,0, //muonID
					     5.,0,          // ET, L1isolation
					     999., 999.,       // Track iso barrel, Track iso endcap
					     999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
					     999., 999.,       // H/ET iso barrel, H/ET iso endcap
					     999., 999.,       // E/ET iso barrel, E/ET iso endcap
					     999., 999.,       // H/E barrel, H/E endcap
					     999., 999.,       // cluster shape barrel, cluster shape endcap
					     0.98, 1.0,       // R9 barrel, R9 endcap
					     999., 999.,       // Deta barrel, Deta endcap
					     999., 999.        // Dphi barrel, Dphi endcap
					     )>1)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }

  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_HT250_DoubleLooseIsoPFTau10trk3_PFMHT35") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(35.) && 
	   OpenHltSumCorHTPassed(250.,40.,3.)>0 && 
	   OpenHltPFTauPassedNoMuonIDNoEleID(10.,3.,1.5,999.,
					     3.,3.,5.,2.,0, //muonID
					     5.,0,          // ET, L1isolation
					     999., 999.,       // Track iso barrel, Track iso endcap
					     999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
					     999., 999.,       // H/ET iso barrel, H/ET iso endcap
					     999., 999.,       // E/ET iso barrel, E/ET iso endcap
					     999., 999.,       // H/E barrel, H/E endcap
					     999., 999.,       // cluster shape barrel, cluster shape endcap
					     0.98, 1.0,       // R9 barrel, R9 endcap
					     999., 999.,       // Deta barrel, Deta endcap
					     999., 999.        // Dphi barrel, Dphi endcap
					     )>1)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }
  
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_HT300_DoubleLooseIsoPFTau10trk3_PFMHT40") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(40.) && 
	   OpenHltSumCorHTPassed(300.,40.,3.)>0 && 
	   OpenHltPFTauPassedNoMuonIDNoEleID(10.,3.,1.5,999.,
					     3.,3.,5.,2.,0, //muonID
					     5.,0,          // ET, L1isolation
					     999., 999.,       // Track iso barrel, Track iso endcap
					     999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
					     999., 999.,       // H/ET iso barrel, H/ET iso endcap
					     999., 999.,       // E/ET iso barrel, E/ET iso endcap
					     999., 999.,       // H/E barrel, H/E endcap
					     999., 999.,       // cluster shape barrel, cluster shape endcap
					     0.98, 1.0,       // R9 barrel, R9 endcap
					     999., 999.,       // Deta barrel, Deta endcap
					     999., 999.        // Dphi barrel, Dphi endcap
					     )>1)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }
  
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_HT350_DoubleLooseIsoPFTau10trk3_PFMHT45") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1){
      if (prescaleResponse(menu,cfg,rcounter,it)){
	if(OpenHltpfMHT(45.) && 
	   OpenHltSumCorHTPassed(350.,40.,3.)>0 && 
	   OpenHltPFTauPassedNoMuonIDNoEleID(10.,3.,1.5,999.,
					     3.,3.,5.,2.,0, //muonID
					     5.,0,          // ET, L1isolation
					     999., 999.,       // Track iso barrel, Track iso endcap
					     999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
					     999., 999.,       // H/ET iso barrel, H/ET iso endcap
					     999., 999.,       // E/ET iso barrel, E/ET iso endcap
					     999., 999.,       // H/E barrel, H/E endcap
					     999., 999.,       // cluster shape barrel, cluster shape endcap
					     0.98, 1.0,       // R9 barrel, R9 endcap
					     999., 999.,       // Deta barrel, Deta endcap
					     999., 999.        // Dphi barrel, Dphi endcap
					     )>1)
	  {
	    triggerBit[it] = true;
	  }
      }
    }
  }

	
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT50U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && 
		  OpenHltSumHTPassed(
				     50,
				     20)>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT70U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && OpenHltSumHTPassed(
										 70,
										 20)>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_MET20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && recoMetCal>45)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_MET45") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && recoMetCal>45)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet30U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1
		  && OpenHlt1JetPassed(30)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet35U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1
		  && OpenHlt1JetPassed(35)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet50U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1
		  && OpenHlt1JetPassed(50)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu8_Jet40_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1JetPassed(40.)>=1)
		{
		  if (OpenHlt1MuonPassed(3., 4., 8., 2., 0)>=1)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet70") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1JetPassed(70.)>=1)
		{
		  if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1)
		    triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_MET45x") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && recoMetCal>45)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu7_MET20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 7., 2., 0)>=1 && recoMetCal>20)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT70U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && OpenHltSumHTPassed(
										 70,
										 20)>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT100U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && OpenHltSumHTPassed(
										 100,
										 20)>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT120U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && OpenHltSumHTPassed(
									     120,
									     20)>0)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_HT140U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1 && OpenHltSumHTPassed(
										 140,
										 20)>0)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu20_CentralJet20U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 7., 20., 2., 0)>=1 && OpenHlt1JetPassed(
										 20,
										 2.6)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu17_TripleCentralJet20U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 7., 17., 2., 0)>=1 && OpenHlt1JetPassed(
										 20,
										 2.6)>=3)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet50U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1
		  && OpenHlt1JetPassed(50) >=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_Jet70U") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=1
		  && OpenHlt1JetPassed(70) >=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
 
	
  else if (triggerName.CompareTo("OpenHLT_L2Mu8_HT50") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(50., 30.) == 1)
		{
		  if (OpenHlt1L2MuonPassed(7., 8., 9999.)>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_L2Mu10_HT50") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second>0)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(50., 30.) == 1)
		{
		  if (OpenHlt1L2MuonPassed(7., 10., 9999.)>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
  /*muon-Tau cross-triggers*/
  else if (triggerName.CompareTo("OpenHLT_Mu11_PFIsoTau15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 7., 11., 2., 0)>=1)
		if (OpenHltPFTauPassedNoMuon(15., 1., 1, 1.)>=1)
		  {
		    triggerBit[it] = true;
		  }
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu15_PFIsoTau15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 7., 15., 2., 0)>=1)
		if (OpenHltPFTauPassedNoMuon(15., 1., 1, 1)>=1)
		  {
		    triggerBit[it] = true;
		  }
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu11_PFIsoTau15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(7., 7., 11., 2., 1)>=1)
		if (OpenHltPFTauPassedNoMuon(15., 1., 1, 1)>=1)
		  {
		    triggerBit[it] = true;
		  }
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu12_PFIsoTau10_Trk1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 12., 2., 1)>=1)
		{
		  if (OpenHltPFTauPassedNoMuon(10., 1., 1., 1.) >=1)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu17_PFIsoTau15_Trk5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(14., 14., 17., 2., 0)>=1)
		{
		  if (OpenHltPFTauPassedNoMuon(15., 5., 1., 1.) >=1)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu15_PFIsoTau20_Trk5") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 15., 2., 1)>=1)
		{
		  if (OpenHltPFTauPassedNoMuon(20., 5., 1., 1.) >=1)
		    triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu15_LooseIsoPFTau15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 15., 2., 0)>=1)
		if (OpenHltPFTauPassedNoMuon(15., 1., 1.5, 999.)>=1)
		  {
		    triggerBit[it] = true;
		  }
	    }
	}
    }

  else if (triggerName.CompareTo("OpenHLT_IsoMu15_LooseIsoPFTau15") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(10., 10., 15., 2., 1)>=1)
		if (OpenHltPFTauPassedNoMuon(15., 1., 1.5, 999.)>=1)
		  {
		    triggerBit[it] = true;
		  }
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu15_LooseIsoPFTau20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(10., 10., 15., 2., 1)>=1)
                if (OpenHltPFTauPassedNoMuon(20., 1., 1.5, 999.)>=1)
                  {
                    triggerBit[it] = true;
                  }
            }
        }
    }
  else if (triggerName.CompareTo("OpenHLT_IsoMu15_MediumIsoPFTau20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1MuonPassed(10., 10., 15., 2., 1)>=1)
                if (OpenHltPFTauPassedNoMuon(20., 1., 1., 999.)>=1)
                  {
                    triggerBit[it] = true;
                  }
            }
        }
    }

  /*Electron-Tau cross-triggers*/
  else if (triggerName.CompareTo("OpenHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          18.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
	      triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          18.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }
	
  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta2p1_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  2.1,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta1p5_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  1.5,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta1p0_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  1.0,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20Eta2p1")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.1,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20Eta1p5")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          1.5,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20Eta1p0")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          1.0,
                                          1.,
                                          1.,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }


  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2")
	   == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1Ele1PFTauPassed(
					  15.,
					  999.,
					  0,
					  999.,
					  999.,
					  999.,
					  999.,
					  999.,
					  999.,
					  999.,
					  999.,
					  0.05,
					  0.05,
					  0.011,
					  0.031,
					  0.98,
					  1.,
					  0.008,
					  0.008,
					  0.07,
					  0.05,
					  15.,
					  2.5,
					  1.,
					  1.5,
					  1000.,
					  0.) >= 1)
		triggerBit[it] = true;
	    }
	}
    }

  else if (triggerName.CompareTo("OpenHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          15.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.5,
                                          1000.,
                                          0.) >= 1)
                triggerBit[it] = true;
            }
        }
    }

	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1Ele1PFTauPassed(
					  15.,
					  999.,
					  0,
					  999.,
					  999.,
					  0.125,
					  0.075,
					  0.125,
					  0.075,
					  0.125,
					  0.075,
					  0.05,
					  0.05,
					  0.011,
					  0.031,
					  0.98,
					  1.,
					  0.008,
					  0.008,
					  0.07,
					  0.05,
					  15.,
					  2.5,
					  1.,
					  1.5,
					  1000.,
					  0.) >= 1)
		triggerBit[it] = true;
	    }
	}
    }
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1Ele1PFTauPassed(
					  15.,
					  999.,
					  0,
					  999.,
					  999.,
					  0.125,
					  0.075,
					  0.125,
					  0.075,
					  0.125,
					  0.075,
					  0.05,
					  0.05,
					  0.011,
					  0.031,
					  0.98,
					  1.,
					  0.008,
					  0.008,
					  0.07,
					  0.05,
					  20.,
					  2.5,
					  1.,
					  1.5,
					  1000.,
					  0.) >= 1)
		triggerBit[it] = true;
	    }
	}
    }

  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          18.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.5,
                                          1000.,
                                          0.) >= 1)
                triggerBit[it] = true;
            }
        }
    }
	
  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
                                          1.5,
                                          1000.,
                                          0.) >= 1)
                triggerBit[it] = true;
            }
        }
    }
	
  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_TrkIdT_TrkIsoT_LooseIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          999.,
                                          999.,
                                          999.,
                                          999.,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta2p1_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  2.1,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta1p5_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  1.5,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20Eta1p0_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  1.0,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.5,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20Eta2p1")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          2.1,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20Eta1p5")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          1.5,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

  else if (triggerName.CompareTo("OpenHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20Eta1p0")
           == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              if (OpenHlt1Ele1PFTauPassed(
                                          20.,
					  999.,
                                          0,
                                          999.,
                                          999.,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.125,
                                          0.075,
                                          0.05,
                                          0.05,
                                          0.011,
                                          0.031,
                                          0.98,
                                          1.,
                                          0.008,
                                          0.008,
                                          0.07,
                                          0.05,
                                          20.,
                                          1.0,
                                          1.,
					  1.5,
                                          1000.,
                                          0.) >= 1)
              triggerBit[it] = true;
            }
        }
    }

	// hartl 2011-06-12, requested by Alexander Gude
	
   else if (triggerName.CompareTo("OpenHLT_Ele17_CaloIdIso_HFEM15_L1R_v1") == 0) {
       if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1) {
           if (prescaleResponse(menu, cfg, rcounter, it)) {
               // 17 GeV Electron
               if ( 
                       OpenHlt1ElectronSamHarperPassed(
                           17., 0, // ET, L1isolation  
                           999., 999., // Track iso barrel, Track iso endcap  
                           999., 999., // Track/pT iso barrel, Track/pT iso endcap  
                           .2, .2, // H/ET iso barrel, H/ET iso endcap  
                           .2, .2, // E/ET iso barrel, E/ET iso endcap  
                           0.15, 0.1,  // H/E barrel, H/E endcap  
                           0.014, 0.035, // cluster shape barrel, cluster shape endcap  
                           0.98, 1.0,  // R9 barrel, R9 endcap  
                           999, 999,   // Deta barrel, Deta endcap  
                           999, 999    // Dphi barrel, Dphi endcap  
                           ) >=1 
                  ) {
                   // HF e code
                   float MinPt  = 15.0;
                   float MaxEta = 5.0;
                   float MinEta = -5.0;
                   float e9e25Min = 0.90; // e9e25 > 0.90 
                   float Min2DCut = 0.2; // 2DCut > 0.2

                   for(int i = 0; i < NohHFECALClus; i++){ // These arrays are kept in sync when creating openhlt.root
                       if (
                               ohHFEleClustere9e25[i] > e9e25Min &&
                               ohHFEleCluster2Dcut[i] > Min2DCut &&
                               ohHFEleEta[i] < MaxEta &&
                               ohHFEleEta[i] > MinEta &&
                               ohHFElePt[i]  > MinPt
                          ) {
                           triggerBit[it] = true; // Fired
                       } 
                   }
               }
           } 
       }
   }

  /* Tau-jet/MET cross-triggers */
  else if (isQuadJetX_IsoPFTauXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltQuadJetPassedPlusTauPFId(
						  thresholds[0],
						  2.5,
						  thresholds[1]) == 1 && OpenL1QuadJet8(10, 2.5) >= 4)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isQuadJetX_IsoPFTauX_PFMHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltQuadJetPassedPlusTauPFId(
						  thresholds[0],
						  2.5,
						  thresholds[1]) == 1 && OpenL1QuadJet8(10, 2.5) >= 4)
		{
		  if (pfMHT > thresholds[2])
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
  else if (menu->GetTriggerName(it).CompareTo("OpenHLT_QuadJet30_NewIsoPFTauTag") == 0) {
    if (map_L1BitOfStandardHLTPath.find(menu->GetTriggerName(it))->second==1) {
      if(OpenHltQuadJetCORPassedPlusTauPFIdNewIso(30, 2.5, 30) == 1 && OpenL1QuadJet8(10, 2.5) >= 4) {
	if (prescaleResponse(menu,cfg,rcounter,it)) { triggerBit[it] = true; }
      }
    }
  }
	
	
  /***********OpenHLT_SingleIsoTauX_TrkX_METX***********/
  else if (isSingleIsoTauX_TrkX_METXTrigger(
					    triggerName,
					    thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltTauL2SCMETPassed(
					  thresholds[0],
					  thresholds[1],
					  0,
					  0.,
					  1,
					  thresholds[2],
					  20.,
					  30.)>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*****************************************************/
	
  else if (triggerName.CompareTo("OpenHLT_SingleTau5_Trk0_MET50_Level1_10") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (recoMetCal>50)
		{
		  //                if(OpenHltTauL2SCMETPassed(5.,0.,0,0.,0,50.,10.,10.)>=1) {
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_SinglePFIsoTau50_Trk15_PFMHT40") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPFTauPassedNoMuon(50., 15., 1, 1)>=1)
		{
		  if (pfMHT > 40)
		    triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_SinglePFIsoTau30_Trk15_PFMHT50") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPFTauPassedNoMuon(30., 15., 1, 1)>=1)
		{
		  if (pfMHT > 50)
		    triggerBit[it] = true;
		}
	    }
	}
    }
	
  /* Electron-MET cross-triggers */
	
  // 2011-03-29 promoted to v3 TODO check
  else if (isEleX_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ( (OpenHlt1ElectronSamHarperPassed(thresholds[0], 0, // ET, L1isolation 
						    999.,
						    999., // Track iso barrel, Track iso endcap 
						    0.2,
						    0.2, // Track/pT iso barrel, Track/pT iso endcap 
						    0.2,
						    0.2, // H/ET iso barrel, H/ET iso endcap 
						    0.2,
						    0.2, // E/ET iso barrel, E/ET iso endcap 
						    0.15,
						    0.10, // H/E barrel, H/E endcap 
						    0.014,
						    0.035, // cluster shape barrel, cluster shape endcap 
						    0.98,
						    1.0, // R9 barrel, R9 endcap 
						    0.01,
						    0.01, // Deta barrel, Deta endcap 
						    0.15,
						    0.10 // Dphi barrel, Dphi endcap 
						    )>=1)
		   && 
		   (OpenHltSumCorHTPassed(thresholds[1], 40)>=1))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v3 TODO check
  else if (isEleX_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if ( (OpenHlt1ElectronSamHarperPassed(thresholds[0], 0, // ET, L1isolation 
						    999.,
						    999., // Track iso barrel, Track iso endcap 
						    0.2,
						    0.2, // Track/pT iso barrel, Track/pT iso endcap 
						    0.2,
						    0.2, // H/ET iso barrel, H/ET iso endcap 
						    0.2,
						    0.2, // E/ET iso barrel, E/ET iso endcap 
						    0.1,
						    0.075, // H/E barrel, H/E endcap 
						    0.011,
						    0.031, // cluster shape barrel, cluster shape endcap 
						    0.98,
						    1.0, // R9 barrel, R9 endcap 
						    0.008,
						    0.008, // Deta barrel, Deta endcap 
						    0.07,
						    0.05 // Dphi barrel, Dphi endcap 
						    )>=1)
		   && 
		   (OpenHltSumCorHTPassed(thresholds[1], 40)>=1))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /* Jet-MET/HT cross-triggers*/
	
  /****CentralJetXU_METX*****/
  else if (isCentralJetXU_METXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1JetPassed(thresholds[0], 2.6)>=1 && recoMetCal
		  >=thresholds[1])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /****CentralJetX_METX*****/
  else if (isCentralJetX_METXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1CorJetPassed(thresholds[0], 2.6)>=1 && recoMetCal
		  >=thresholds[1])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  /********DiJetX(U)_METX***********************/
  else if (isDiJetXU_METXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1JetPassed(thresholds[0])>=2 && recoMetCal
		  >=thresholds[1])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDiJetX_METXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1CorJetPassed(thresholds[0])>=2 && recoMetCal
		  >=thresholds[1])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  /**DiJetX(U)_PTX(U)**/
	
  else if (isDiJetXU_PTXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPT12U(thresholds[1], thresholds[0]) == 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isDiJetX_PTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltPT12(thresholds[1], thresholds[0]) == 1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /***********/
	
  /*HT-MET/MHT cross-triggers*/
  /*****METX_HTX(U)**************/
  else if (isMETX_HTXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumHTPassed(thresholds[1], 20)>=1 && recoMetCal
		  >=thresholds[0])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isMETX_HTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(thresholds[1], 40)>=1 && recoMetCal
		  >=thresholds[0])
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /***HTX(U)_MHT(X)U*****/
	
  else if (isHTXU_MHTXUTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltMHTU(thresholds[1], 20.)==1 && (OpenHltSumHTPassed(
									    thresholds[0],
									    20.) == 1))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isHTX_MHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltMHT(thresholds[1], 30.)==1 && (OpenHltSumCorHTPassed(
									      thresholds[0],
									      40.) == 1))
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

 else if (isDiCentralPFJetX_PFMHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      if  ( OpenHltNPFJetPassed(2, thresholds[0], 2.6) && OpenHltPFMHT(thresholds[1], 0.)==1)
		  
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }



 else if (isPFHTX_PFMHTX_OrHTX_OrMHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      if (
		  (OpenHltMHT(thresholds[1], 30.)==1 && OpenHltSumCorHTPassed(thresholds[0], 40.) == 1)
		  &&( 
		     (OpenHltPFMHT(thresholds[1], 30.)==1 && OpenHltSumPFHTPassed(thresholds[0], 40.) == 1)
		     || 
		     (OpenHltSumCorHTPassed(thresholds[2], 40.) && OpenHltMHT(thresholds[3], 30.)==1  )
		     )
		  )
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }

  else if (isPFHTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {

	      if (
		  (OpenHltSumCorHTPassed(thresholds[0], 40.) == 1)
		  && OpenHltSumPFHTPassed(thresholds[0], 40.) == 1)
		   
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
	
  /*****HTX_AlphaT0pX*******/
  else if (isAlphaTTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltHT_AlphaT( thresholds[0],thresholds[1],40. ) >=1 )
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
	
  /*********************/
	
  /*Muon-photon cross-triggers*/
	
  else if (triggerName.CompareTo("OpenHLT_Mu8_Photon20_CaloIdVT_IsoT_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 4., 8., 2., 0)>=1
		  && OpenHlt1PhotonSamHarperPassed(20., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   3.0,
						   3.0, // Track/pT iso barrel, Track/pT iso endcap
						   3.0,
						   3.0, // H/ET iso barrel, H/ET iso endcap
						   5.0,
						   5.0, // E/ET iso barrel, E/ET iso endcap
						   0.05,
						   0.05, // H/E barrel, H/E endcap  
						   0.011,
						   0.031, // cluster shape barrel, cluster shape endcap  
						   0.98,
						   999., // R9 barrel, R9 endcap
						   999.,
						   999., // Deta barrel, Deta endcap
						   999.,
						   999. // Dphi barrel, Dphi endcap
						   )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu15_Photon20_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(5., 7., 15., 2., 0)>=1
		  && OpenHlt1PhotonSamHarperPassed(20., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   999.,
						   999., // Track/pT iso barrel, Track/pT iso endcap
						   999.,
						   999., // H/ET iso barrel, H/ET iso endcap
						   999.,
						   999., // E/ET iso barrel, E/ET iso endcap
						   0.15,
						   0.10, // H/E barrel, H/E endcap
						   0.014,
						   0.035, // cluster shape barrel, cluster shape endcap
						   0.98,
						   999., // R9 barrel, R9 endcap
						   999.,
						   999., // Deta barrel, Deta endcap
						   999.,
						   999. // Dphi barrel, Dphi endcap
						   )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  //else if (triggerName.CompareTo("OpenHLT_Mu15_DiPhoton15_CaloIdL_v1") == 0)
  else if (triggerName.CompareTo("OpenHLT_Mu15_DoublePhoton15_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(5., 7., 15., 2., 0)>=1
		  && OpenHlt1PhotonSamHarperPassed(15., 0, // ET, L1isolation
						   999.,
						   999., // Track iso barrel, Track iso endcap
						   999.,
						   999., // Track/pT iso barrel, Track/pT iso endcap
						   999.,
						   999., // H/ET iso barrel, H/ET iso endcap
						   999.,
						   999., // E/ET iso barrel, E/ET iso endcap
						   0.15,
						   0.10, // H/E barrel, H/E endcap
						   0.014,
						   0.035, // cluster shape barrel, cluster shape endcap
						   0.98,
						   999., // R9 barrel, R9 endcap
						   999.,
						   999., // Deta barrel, Deta endcap
						   999.,
						   999. // Dphi barrel, Dphi endcap
						   )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  /*Muon-electron cross-triggers*/
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_DoubleEle8_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=1
		  &&OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.1, // H/E barrel, H/E endcap 
						     999.,
						     999., // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1., // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=2)
		{
		  //OpenHlt1ElectronPassed(8.,0,9999.,9999.)>=2){ 
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=1
		  &&OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.1, // H/E barrel, H/E endcap 
						     999.,
						     999., // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1., // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=2 &&
		  OpenHlt1ElectronSamHarperPassed(8.,0,          // ET, L1isolation
						  999., 999.,       // Track iso barrel, Track iso endcap
						  999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
						  999., 999.,       // H/ET iso barrel, H/ET iso endcap
						  999., 999.,       // E/ET iso barrel, E/ET iso endcap
						  0.15, 0.10,       // H/E barrel, H/E endcap 
						  0.014, 0.035,       // cluster shape barrel, cluster shape endcap 
						  0.98, 1.,       // R9 barrel, R9 endcap
						  0.01, 0.01,       // Deta barrel, Deta endcap 
						  0.15, 0.10        // Dphi barrel, Dphi endcap
						  )>=1)
		{
		  //OpenHlt1ElectronPassed(8.,0,9999.,9999.)>=2){ 
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  else if (triggerName.CompareTo("OpenHLT_DoubleMu5_Ele8_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(0., 0., 5., 2., 0)>=2
		  && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.1, // H/E barrel, H/E endcap 
						     999.,
						     999., // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1.0, // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 4., 5., 2., 0)>=2
		  && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.10, // H/E barrel, H/E endcap  
						     0.014,
						     0.035, // cluster shape barrel, cluster shape endcap  
						     0.98,
						     1., // R9 barrel, R9 endcap  
						     0.01,
						     0.01, // Deta barrel, Deta endcap  
						     0.15,
						     0.10 // Dphi barrel, Dphi endcap 
						     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (isMuX_EleX_HTXTrigger(triggerName, thresholds, caloId, caloIso, trkId, trkIso)){
    
    if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
	if (prescaleResponse(menu, cfg, rcounter, it))
	  {
	    if ((OpenHlt1MuonPassed(map_muThresholds[thresholds[0]], 2., 0)>=1) &&
	        (OpenHlt1ElectronPassed(thresholds[1], 
				       map_EGammaCaloId[caloId[0]],
				       map_EleCaloIso[caloIso[0]],
				       map_EleTrkId[trkId[0]],
				       map_EleTrkIso[trkIso[0]]
					) >= 1) &&
		(OpenHltSumCorHTPassed(thresholds[2], 40.) == 1))
	      triggerBit[it] = true;
	      }
	  }
      }


  else if (triggerName.CompareTo("OpenHLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 3., 2., 0)>=1 && OpenHltSumCorHTPassed(
										    160,
										    40)>0 && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation
															     999.,
															     999., // Track iso barrel, Track iso endcap
															     999.,
															     999., // Track/pT iso barrel, Track/pT iso endcap
															     999.,
															     999., // H/ET iso barrel, H/ET iso endcap
															     999.,
															     999., // E/ET iso barrel, E/ET iso endcap
															     0.15,
															     0.10, // H/E barrel, H/E endcap
															     0.014,
															     0.035, // cluster shape barrel, cluster shape endcap
															     0.98,
															     1.0, // R9 barrel, R9 endcap
															     0.01,
															     0.01, // Deta barrel, Deta endcap
															     0.15,
															     0.10 // Dphi barrel, Dphi endcap
															     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
  else if (triggerName.CompareTo("OpenHLT_Mu5_DoubleEle8_CaloIdL_TrkIdVL") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=1
		  &&OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15, 0.10,       // H/E barrel, H/E endcap 
						     0.014, 0.035,       // cluster shape barrel, cluster shape endcap 
						     0.98, 1.,       // R9 barrel, R9 endcap
						     0.01, 0.01,       // Deta barrel, Deta endcap 
						     0.15, 0.10        // Dphi barrel, Dphi endcap
						     )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_CaloIdT_TrkIdVL for Jonathan
  else if (triggerName.CompareTo("OpenHLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_CaloIdT_TrkIdVL") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 5., 2., 0)>=1
		  &&OpenHlt2ElectronsSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15, 0.10,       // H/E barrel, H/E endcap 
						     0.014, 0.035,       // cluster shape barrel, cluster shape endcap 
						     0.98, 1.,       // R9 barrel, R9 endcap
						     0.01, 0.01,       // Deta barrel, Deta endcap 
						     0.15, 0.10        // Dphi barrel, Dphi endcap
						     )>=2 &&
		  OpenHlt1ElectronSamHarperPassed(8.,0,          // ET, L1isolation
						  999., 999.,       // Track iso barrel, Track iso endcap
						  999., 999.,        // Track/pT iso barrel, Track/pT iso endcap
						  999., 999.,       // H/ET iso barrel, H/ET iso endcap
						  999., 999.,       // E/ET iso barrel, E/ET iso endcap
						  0.1, 0.075, // H/E barrel, H/E endcap
						  0.011, 0.031, // cluster shape barrel, cluster shape endcap
						  0.98, 1.,       // R9 barrel, R9 endcap
						  0.01, 0.01,       // Deta barrel, Deta endcap 
						  0.15, 0.10        // Dphi barrel, Dphi endcap
						  )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
	
  else if (triggerName.CompareTo("OpenHLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 3., 3., 2., 0)>=1 && OpenHltSumCorHTPassed(
										    160,
										    40)>0 && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation
															     999.,
															     999., // Track iso barrel, Track iso endcap
															     999.,
															     999., // Track/pT iso barrel, Track/pT iso endcap
															     999.,
															     999., // H/ET iso barrel, H/ET iso endcap
															     999.,
															     999., // E/ET iso barrel, E/ET iso endcap
															     0.10,
															     0.075, // H/E barrel, H/E endcap 
															     0.011,
															     0.031, // cluster shape barrel, cluster shape endcap 
															     0.98,
															     1.0, // R9 barrel, R9 endcap
															     0.01,
															     0.01, // Deta barrel, Deta endcap
															     0.15,
															     0.10 // Dphi barrel, Dphi endcap
															     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v3 TODO check
  else if (isDoubleEleX_CaloIdL_TrkIdVL_HTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(thresholds[1], 40)>0
		  && OpenHlt2ElectronsSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						      999.,
						      999., // Track iso barrel, Track iso endcap
						      999.,
						      999., // Track/pT iso barrel, Track/pT iso endcap
						      999.,
						      999., // H/ET iso barrel, H/ET iso endcap
						      999.,
						      999., // E/ET iso barrel, E/ET iso endcap
						      0.15,
						      0.1, // H/E barrel, H/E endcap
						      0.014,
						      0.035, // cluster shape barrel, cluster shape endcap
						      0.98,
						      1.0, // R9 barrel, R9 endcap
						      0.1,
						      0.1, // Deta barrel, Deta endcap
						      0.15,
						      0.1 // Dphi barrel, Dphi endcap
						      )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  // 2011-03-29 promoted to v3 TODO check
  else if (isDoubleEleX_CaloIdT_TrkIdVL_HTXTrigger(triggerName, thresholds))
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHltSumCorHTPassed(thresholds[1], 40)>0
		  && OpenHlt2ElectronsSamHarperPassed(thresholds[0], 0, // ET, L1isolation
						      999.,
						      999., // Track iso barrel, Track iso endcap
						      999.,
						      999., // Track/pT iso barrel, Track/pT iso endcap
						      999.,
						      999., // H/ET iso barrel, H/ET iso endcap
						      999.,
						      999., // E/ET iso barrel, E/ET iso endcap
						      0.1,
						      0.075, // H/E barrel, H/E endcap
						      0.011,
						      0.031, // cluster shape barrel, cluster shape endcap
						      0.98,
						      1.0, // R9 barrel, R9 endcap
						      0.1,
						      0.1, // Deta barrel, Deta endcap
						      0.15,
						      0.1 // Dphi barrel, Dphi endcap
						      )>=2)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu17_Ele8_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 5., 17., 2., 0)>=1
		  && OpenHlt1ElectronSamHarperPassed(8., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.10, // H/E barrel, H/E endcap 
						     0.014,
						     0.035, // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1.0, // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu10_Ele10_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 5., 10., 2., 0)>=1
		  && OpenHlt1ElectronSamHarperPassed(10., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.10, // H/E barrel, H/E endcap 
						     0.014,
						     0.035, // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1.0, // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	
  else if (triggerName.CompareTo("OpenHLT_Mu8_Ele17_CaloIdL_v1") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      if (OpenHlt1MuonPassed(3., 5., 8., 2., 0)>=1
		  && OpenHlt1ElectronSamHarperPassed(17., 0, // ET, L1isolation 
						     999.,
						     999., // Track iso barrel, Track iso endcap 
						     999.,
						     999., // Track/pT iso barrel, Track/pT iso endcap 
						     999.,
						     999., // H/ET iso barrel, H/ET iso endcap 
						     999.,
						     999., // E/ET iso barrel, E/ET iso endcap 
						     0.15,
						     0.10, // H/E barrel, H/E endcap 
						     0.014,
						     0.035, // cluster shape barrel, cluster shape endcap 
						     0.98,
						     1.0, // R9 barrel, R9 endcap 
						     999.,
						     999., // Deta barrel, Deta endcap 
						     999.,
						     999. // Dphi barrel, Dphi endcap 
						     )>=1)
		{
		  triggerBit[it] = true;
		}
	    }
	}
    }
	 
	 
	
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon36_IsoVL_Photon22_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									22.,  
									0,
									999., //was 4
									999., //was 6
									999.,  //was 4
									999.,  //was 4
									0.15,   
									0.10,   
									999.); //was 0.98
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
									     36.,  
									     0,
									     4.0,   
									     6.0,   
									     4.0,   
									     4.0,   
									     0.15,  
									     0.10,  
									     999.); 
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  // 2011-04-26
  else if (triggerName.CompareTo("OpenHLT_Photon40_CaloIdL_Photon28_CaloIdL") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									28.,
									0,
									999,
									999,
									999,
									999,
									0.15,
									0.10,
									999.,
									0.014,
									0.035);
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
									     40., 
									     0,
									     999,
									     999,
									     999,
									     999,
									     0.15, 
									     0.10, 
									     999.,
									     0.014, 
									     0.035); 
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  // 2011-04-26 
  else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdL_IsoVL_Photon22_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									22., 
									0,
									999.,
									999.0,
									999.0,
									999.0,
									0.15, 
									0.10, 
									999.);
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
									     36., 
									     0,
									     4.0, 
									     6.0, 
									     4.0, //HBR  
									     4.0, //HEC  
									     0.15,  
									     0.10,  
									     999.,
									     0.014, 
									     0.035); 
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  // 2011-04-26 
  else if (triggerName.CompareTo("OpenHLT_HLT_Photon26_IsoL_Photon18_v3") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									18., 
									0,
									999.,
									999.0,
									999.0,
									999.0,
									0.15, 
									0.10, 
									999.); 
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
									     26., 
									     0,
									     3.5, 
									     5.5,  
									     3.5, //HBR   
									     3.5, //HEC   
									     0.15,   
									     0.10,  
									     999.,
									     0.014, 
									     0.035); 
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  // 2011-03-29 promoted to v2 TODO check
  else if (triggerName.CompareTo("OpenHLT_Photon20_R9Id_Photon18_R9Id_v2") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
									    18,
									    0.8,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15, //H/E EB
									    0.10, // H/E EC
									    999.);
				
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector =
		    VectorOpenHlt1PhotonPassedR9ID(
						   20,
						   0.8,
						   0,
						   999,
						   999,
						   999,
						   999,
						   0.15, //H/E EB
						   0.10, //H/E EC
						   999.);
					
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
//2011-06-08
  else if (triggerName.CompareTo("OpenHLT_Photon26_R9Id_Photon18_R9Id_v3") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
        {
          if (prescaleResponse(menu, cfg, rcounter, it))
            {
              std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
                                                                            18,
                                                                            0.8,
                                                                            0,
                                                                            999,
                                                                            999,
                                                                            999,
                                                                            999,
                                                                            0.15, //H/E EB
                                                                            0.10, // H/E EC
                                                                            999.);

              if (firstVector.size()>=2)
                {
                  std::vector<int> secondVector =
                    VectorOpenHlt1PhotonPassedR9ID(
                                                   26,
                                                   0.8,
                                                   0,
                                                   999,
                                                   999,
                                                   999,
                                                   999,
                                                   0.15, //H/E EB
                                                   0.10, //H/E EC
                                                   999.);

                  if (secondVector.size()>=1)
                    {
                      triggerBit[it] = true;
                    }
                }
            }
        }
    }
	
  // 2011-05-04
  else if (triggerName.CompareTo("OpenHLT_Photon26_R9ID_Photon18_R9ID_v1") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
									    18,   
									    0.8,  
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15, //H/E EB     
									    0.10, // H/E EC
									    999.);              
				
	      if (firstVector.size()>=2)
		{
		  std::vector<int> secondVector =
		    VectorOpenHlt1PhotonPassedR9ID(
						   26,  
						   0.8,  
						   0,
						   999,
						   999,
						   999,
						   999,
						   0.15, //H/E EB
						   0.10, //H/E EC   
						   999.);            
					
		  if (secondVector.size()>=1)
		    {
		      triggerBit[it] = true;
		    }
		}
	    }
	}
    }
	
	
  // 2011-05-04
  else if (triggerName.CompareTo("OpenHLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
									    18,//Et  
									    0.8, //R9ID 
									    0,
									    999.0, //TrkIso 
									    6.0, //ECAL Iso 
									    4.0,//HCAL Iso BR 
									    4.0, //HCAL Iso EC 
									    0.15, //H/E EB 
									    0.10, // H/E EC 
									    999.,//R9
									    0.014,//Cluster Shape EB 
									    0.035);//Cluster Shape EC
				
	      if (firstVector.size()>=1)
		{
		  std::vector<int> secondVector = VectorOpenHlt1PhotonPassedR9ID(
										 26,//Et 
										 0.,//R9ID 
										 0,
										 4.0,//TrkIso 
										 6.0,//ECALIso 
										 4.0,//HcalIso Br 
										 4.0,//HcalIso EC
										 0.15, //H/E EB 
										 0.10, // H/E EC
										 999.,
										 0.014, //Cluster Shape EB 
										 0.035); //Cluster Shape EC 
					
		  if (secondVector.size()>=1)
		    {
		      bool foundone=false;
		      for (unsigned int j=0; j<firstVector.size() && !foundone; j++) {
			for (unsigned int k=0; k<secondVector.size(); k++) {
			  if (firstVector[j]!=secondVector[k]) {  //if it's not the same object in the ohPhot list
			    //compute dR(j,k)
			    float deltaEta = ohPhotEta[firstVector[j]]-ohPhotEta[secondVector[k]];
			    float deltaPhi = TMath::Abs(ohPhotPhi[firstVector[j]]-ohPhotPhi[secondVector[k]]);
			    if (deltaPhi>TMath::Pi())
			      deltaPhi = 2.0*TMath::Pi()-deltaPhi;
			    double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
									
			    if (deltaR>0.01) {
			      foundone=true;
			      triggerBit[it] = true;
			      break;
			    }
			  }
			}//end for k
		      }//end for j
		    }//end if second vector has any
		}//end if photon 18
	    }
	}
    }
	
	
  // 2011-05-04
  else if (triggerName.CompareTo("OpenHLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3") == 0)//new
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
									    18,   
									    0.8,
									    0,
									    4.0,
									    6.0,
									    4.0,
									    4.0,
									    0.15, //H/E EB     
									    0.10, // H/E EC    
									    999.,
									    0.014, 
									    0.035); 
	      if (firstVector.size()>=1)
		{
		  std::vector<int> secondVector =
		    VectorOpenHlt1PhotonPassedR9ID(
						   26,  
						   0.8,  
						   0,
						   999.,
						   999.0,
						   999.0,
						   999.0,
						   0.15, //H/E EB 
						   0.10, //H/E EC   
						   999.);
					
		  if (secondVector.size()>=1)
		    {
		      bool foundone=false;
		      for (unsigned int j=0; j<firstVector.size() && !foundone; j++) {
			for (unsigned int k=0; k<secondVector.size(); k++) {
			  if (firstVector[j]!=secondVector[k]) { //if it's not the same object in the ohPhot list
			    //compute dR(j,k)
			    float deltaEta = ohPhotEta[firstVector[j]]-ohPhotEta[secondVector[k]];
			    float deltaPhi = TMath::Abs(ohPhotPhi[firstVector[j]]-ohPhotPhi[secondVector[k]]);
			    if (deltaPhi>TMath::Pi())
			      deltaPhi = 2.0*TMath::Pi()-deltaPhi;
			    double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
									
			    if (deltaR>0.01) {
			      foundone=true;
			      triggerBit[it] = true;
			      break;
			    }
			  }
			}//end for k
		      }//end for j
		    }//if pho26
		}//if pho18
	    }
	}
    }
 
   // 2011-05-23 hartl: CHECK
   else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
            22,
            0, //L1iso
            999., //Tiso
            999., //Eiso
            999., //HisoBR
            999., //HisoEC
            0.15, //HoverEEB
            0.10, //HoverEEC
            999.,//0.98, //R9
            0.014, //ClusShapEB
            0.035); //ClusShapEC
            if (firstVector.size()>=2)
            {
               std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
               36,
               0, //L1iso
               4.0, //Tiso
               6.0, //Eiso
               4.0, //HisoBR
               4.0, //HisoEC
               0.15, //HoverEEB
               0.10, //HoverEEC
               999.,//0.98, //R9
               0.014, //ClusShapEB
               0.035); //ClusShapEC
               if (secondVector.size()>=1)
               {
                  triggerBit[it] = true;
               }
            }
         }
      }
   }

   // 2011-05-23 hartl
   else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
            22,
            0, //L1iso
            4.0, //Tiso
            6.0, //Eiso
            4.0, //HisoBR
            4.0, //HisoEC
            0.15, //HoverEEB
            0.10, //HoverEEC
            999.,//0.98, //R9
            0.014, //ClusShapEB
            0.035); //ClusShapEC
            if (firstVector.size()>=2)
            {
               std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
               36,
               0, //L1iso
               4.0, //Tiso
               6.0, //Eiso
               4.0, //HisoBR
               4.0, //HisoEC
               0.15, //HoverEEB
               0.10, //HoverEEC
               999.,//0.98, //R9
               0.014, //ClusShapEB
               0.035); //ClusShapEC
               if (secondVector.size()>=1)
               {
                  triggerBit[it] = true;
               }
            }
         }
      }
   }
  
   // 2011-05-23 hartl: CHECK
   else if (triggerName.CompareTo("OpenHLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id") == 0)//new

   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
            22,//Et  
            0.8, //R9ID 
            0,
            999., //TrkIso 
            999., //ECAL Iso 
            999.,//HCAL Iso BR 
            999., //HCAL Iso EC 
            0.15, //H/E EB 
            0.10, // H/E EC 
            999.,//R9
            999.,//Cluster Shape EB 
            999.);//Cluster Shape EC

            if (firstVector.size()>=1)
            {
               std::vector<int> secondVector = VectorOpenHlt1PhotonPassedR9ID(
               36,//Et 
               0.0,//R9ID 
               0,
               4.0,//TrkIso 
               6.0,//ECALIso 
               4.0,//HcalIso Br 
               4.0,//HcalIso EC
               0.15, //H/E EB 
               0.10, // H/E EC
               999.,
               0.014, //Cluster Shape EB 
               0.035); //Cluster Shape EC 

               if (secondVector.size()>=1)
               {
                  bool foundone=false;
                  for (unsigned int j=0; j<firstVector.size() && !foundone; j++)
                  {
                     for (unsigned int k=0; k<secondVector.size(); k++)
                     {
                        if (firstVector[j]!=secondVector[k])
                        { //if it's not the same object in the ohPhot list
                           //compute dR(j,k)
                           float deltaEta = ohPhotEta[firstVector[j]]-ohPhotEta[secondVector[k]];
                           float deltaPhi = TMath::Abs(ohPhotPhi[firstVector[j]]-ohPhotPhi[secondVector[k]]);
                           if (deltaPhi>TMath::Pi())
                           deltaPhi = 2.0*TMath::Pi()-deltaPhi;
                           double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));

                           if (deltaR>0.01)
                           {
                              foundone=true;
                              triggerBit[it] = true;
                              break;
                           }
                        }
                     }//end for k
                  }//end for j
               }//end if second vector has any
            }//end if photon 22
         }
      }
   }

   // 2011-05-23 hartl: CHECK
   else if (triggerName.CompareTo("OpenHLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL") == 0)//new

   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
            22,
            0.0,
            0,
            4.0,
            6.0,
            4.0,
            4.0,
            0.15, //H/E EB     
            0.10, // H/E EC    
            999.,
            0.014,
            0.035);
            if (firstVector.size()>=1)
            {
               std::vector<int> secondVector =
               VectorOpenHlt1PhotonPassedR9ID(
               36,
               0.8,
               0,
               999.,
               999.0,
               999.0,
               999.0,
               0.15, //H/E EB 
               0.10, //H/E EC   
               999.);

               if (secondVector.size()>=1)
               {
                  bool foundone=false;
                  for (unsigned int j=0; j<firstVector.size() && !foundone; j++)
                  {
                     for (unsigned int k=0; k<secondVector.size(); k++)
                     {
                        if (firstVector[j]!=secondVector[k])
                        { //if it's not the same object in the ohPhot list
                           //compute dR(j,k)
                           float deltaEta = ohPhotEta[firstVector[j]]-ohPhotEta[secondVector[k]];
                           float deltaPhi = TMath::Abs(ohPhotPhi[firstVector[j]]-ohPhotPhi[secondVector[k]]);
                           if (deltaPhi>TMath::Pi())
                           deltaPhi = 2.0*TMath::Pi()-deltaPhi;
                           double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));

                           if (deltaR>0.01)
                           {
                              foundone=true;
                              triggerBit[it] = true;
                              break;
                           }
                        }
                     }//end for k
                  }//end for j
               }//if pho36
            }//if pho22
         }
      }
   }
  
  // 2011-05-23 hartl: CHECK
   else if (triggerName.CompareTo("OpenHLT_Photon36_R9Id_Photon22_R9Id") == 0)
   {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
      {
         if (prescaleResponse(menu, cfg, rcounter, it))
         {
            std::vector<int> firstVector = VectorOpenHlt1PhotonPassedR9ID(
            22,
            0.8,
            0,
            999,
            999,
            999,
            999,
            0.15, //H/E EB
            0.10, // H/E EC
            0.98);

            if (firstVector.size()>=2)
            {
               std::vector<int> secondVector =
               VectorOpenHlt1PhotonPassedR9ID(
               36,
               0.8,
               0,
               999,
               999,
               999,
               999,
               0.15, //H/E EB
               0.10, //H/E EC
               0.98);

               if (secondVector.size()>=1)
               {
                  triggerBit[it] = true;
               }
            }
         }
      }
   }
  
  else if (triggerName.BeginsWith("OpenHLT_Photon26")==1)
    {
      // Photon Paths (V. Rekovic)
		
      int lowerEt = 18;
      int upperEt = 26;
		
      // Names of Photon Paths With Mass Cut
      // -------------------------------------
		
      // No special cuts
      // 2011-03-29 promoted to v2 TODO check
      char pathName_Photon_Photon[100];
      snprintf(
	      pathName_Photon_Photon,
              100,
	      "OpenHLT_Photon%d_Photon%d",
	      upperEt,
	      lowerEt);
		
      // Both legs IsoVL
      // 2011-03-29 promoted to v2 TODO check
      char pathName_Photon_IsoVL_Photon_IsoVL[100];
      snprintf(
	      pathName_Photon_IsoVL_Photon_IsoVL,
              100,
	      "OpenHLT_Photon%d_IsoVL_Photon%d_IsoVL",
	      upperEt,
	      lowerEt);
		
      // One leg IsoVL
      // 2011-03-29 promoted to v2 TODO check
      char pathName_Photon_IsoVL_Photon[100];
      snprintf(
	      pathName_Photon_IsoVL_Photon,
              100,
	      "OpenHLT_Photon%d_IsoVL_Photon%d",
	      upperEt,
	      lowerEt);
		
      // One leg IsoT + Mass>60
      char pathName_Photon_IsoT_Photon_Mass60[100];
      snprintf(
	      pathName_Photon_IsoT_Photon_Mass60,
              100,
	      "OpenHLT_Photon%d_IsoT_Photon%d_Mass60",
	      upperEt,
	      lowerEt);
		
      // Both legs IsoT  + Mass>60
      char pathName_Photon_IsoT_Photon_IsoT_Mass60[100];
      snprintf(
	      pathName_Photon_IsoT_Photon_IsoT_Mass60,
              100,
	      "OpenHLT_Photon%d_IsoT_Photon%d_IsoT_Mass60",
	      upperEt,
	      lowerEt);
		
      // Both legs IsoT
      char pathName_Photon_IsoT_Photon_IsoT[100];
      snprintf(
	      pathName_Photon_IsoT_Photon_IsoT,
              100,
	      "OpenHLT_Photon%d_IsoT_Photon%d_IsoT",
	      upperEt,
	      lowerEt);
		
      // One leg IsoT
      char pathName_Photon_IsoT_Photon[100];
      snprintf(
	      pathName_Photon_IsoT_Photon,
              100,
	      "OpenHLT_Photon%d_IsoT_Photon%d",
	      upperEt,
	      lowerEt);
		
      // One leg IsoL
      char pathName_Photon_IsoL_Photon[100];
      snprintf(
	      pathName_Photon_IsoL_Photon,
              100,
	      "OpenHLT_Photon%d_IsoL_Photon%d",
	      upperEt,
	      lowerEt);
		
      // One leg CaloIdL
      char pathName_Photon_CaloIdL_Photon[100];
      snprintf(
	      pathName_Photon_CaloIdL_Photon,
              100,
	      "OpenHLT_Photon%d_CaloIdL_Photon%d",
	      upperEt,
	      lowerEt);
		
      // Both legs CaloIdL + IsoVL
      // 2011-05-11 promoted to v4 TODO check
      char pathName_Photon_CaloIdL_IsoVL_Photon_CaloIdL_IsoVL[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoVL_Photon_CaloIdL_IsoVL,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoVL_Photon%d_CaloIdL_IsoVL_v4",
	      upperEt,
	      lowerEt);
		
      // One leg CaloIdL + IsoVL
      char pathName_Photon_CaloIdL_IsoVL_Photon[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoVL_Photon,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoVL_Photon%d",
	      upperEt,
	      lowerEt);
		
      // One leg CaloIdL + IsoT + Mass>60
      char pathName_Photon_CaloIdL_IsoT_Photon_Mass60[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoT_Photon_Mass60,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoT_Photon%d_Mass60",
	      upperEt,
	      lowerEt);
		
      // Both legs CaloIdL + IsoT + Mass>60
      char pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT_Mass60[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT_Mass60,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoT_Photon%d_CaloIdL_IsoT_Mass60",
	      upperEt,
	      lowerEt);
		
      // Both legs CaloIdL + IsoT
      char pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoT_Photon%d_CaloIdL_IsoT",
	      upperEt,
	      lowerEt);
		
      // One leg CaloIdL + IsoT
      char pathName_Photon_CaloIdL_IsoT_Photon[100];
      snprintf(
	      pathName_Photon_CaloIdL_IsoT_Photon,
              100,
	      "OpenHLT_Photon%d_CaloIdL_IsoT_Photon%d",
	      upperEt,
	      lowerEt);
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      if (triggerName.CompareTo(pathName_Photon_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 999,
										 999,
										 999,
										 999,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoVL_Photon_IsoVL) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    4.0,
									    6.0,
									    4.0,
									    4.0,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 4.0,
										 6.0,
										 4.0,
										 4.0,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoVL_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 4.0,
										 6.0,
										 4.0,
										 4.0,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoT_Photon_Mass60) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
							
			  TLorentzVector e1;
			  TLorentzVector e2;
			  TLorentzVector meson;
			  float mass = 0.;
			  for (unsigned int i=0; i<firstVector.size(); i++)
			    {
			      for (unsigned int j=0; j<secondVector.size() ; j++)
				{
									
				  if (firstVector[i] == secondVector[j])
				    continue;
				  e1.SetPtEtaPhiM(
						  ohPhotEt[firstVector[i]],
						  ohPhotEta[firstVector[i]],
						  ohPhotPhi[firstVector[i]],
						  0.);
				  e2.SetPtEtaPhiM(
						  ohPhotEt[secondVector[j]],
						  ohPhotEta[secondVector[j]],
						  ohPhotPhi[secondVector[j]],
						  0.);
				  meson = e1 + e2;
				  mass = meson.M();
									
				  if (mass>60)
				    triggerBit[it] = true;
									
				}
			    }
							
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoT_Photon_IsoT_Mass60) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    3,
									    5,
									    3,
									    3,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
							
			  TLorentzVector e1;
			  TLorentzVector e2;
			  TLorentzVector meson;
			  float mass = 0.;
			  for (unsigned int i=0; i<firstVector.size(); i++)
			    {
			      for (unsigned int j=0; j<secondVector.size() ; j++)
				{
									
				  if (firstVector[i] == secondVector[j])
				    continue;
				  e1.SetPtEtaPhiM(
						  ohPhotEt[firstVector[i]],
						  ohPhotEta[firstVector[i]],
						  ohPhotPhi[firstVector[i]],
						  0.);
				  e2.SetPtEtaPhiM(
						  ohPhotEt[secondVector[j]],
						  ohPhotEta[secondVector[j]],
						  ohPhotPhi[secondVector[j]],
						  0.);
				  meson = e1 + e2;
				  mass = meson.M();
									
				  if (mass>60)
				    triggerBit[it] = true;
									
				}
			    }
							
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoT_Photon_IsoT) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    3,
									    5,
									    3,
									    3,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoT_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_IsoL_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3.5,
										 5.5,
										 3.5,
										 3.5,
										 0.15,
										 0.10,
										 0.98);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 999,
										 999,
										 999,
										 999,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
	//2011-05-11 removed R9 cuts.
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoVL_Photon_CaloIdL_IsoVL) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0, //L1iso
									    4.0, //Tiso
									    6.0, //Eiso
									    4.0, //HisoBR
									    4.0, //HisoEC
									    0.15, //HoverEEB
									    0.10, //HoverEEC
									    999.,//0.98, //R9
									    0.014, //ClusShapEB
									    0.035); //ClusShapEC
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0, //L1iso
										 4.0, //Tiso
										 6.0, //Eiso
										 4.0, //HisoBR
										 4.0, //HisoEC
										 0.15, //HoverEEB
										 0.10, //HoverEEC
										 999.,//0.98, //R9
										 0.014, //ClusShapEB
										 0.035); //ClusShapEC
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoVL_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 4.0,
										 6.0,
										 4.0,
										 4.0,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoT_Photon_Mass60) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
							
			  TLorentzVector e1;
			  TLorentzVector e2;
			  TLorentzVector meson;
			  float mass = 0.;
			  for (unsigned int i=0; i<firstVector.size(); i++)
			    {
			      for (unsigned int j=0; j<secondVector.size() ; j++)
				{
									
				  if (firstVector[i] == secondVector[j])
				    continue;
				  e1.SetPtEtaPhiM(
						  ohPhotEt[firstVector[i]],
						  ohPhotEta[firstVector[i]],
						  ohPhotPhi[firstVector[i]],
						  0.);
				  e2.SetPtEtaPhiM(
						  ohPhotEt[secondVector[j]],
						  ohPhotEta[secondVector[j]],
						  ohPhotPhi[secondVector[j]],
						  0.);
				  meson = e1 + e2;
				  mass = meson.M();
									
				  if (mass>60)
				    triggerBit[it] = true;
									
				}
			    }
							
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT_Mass60) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    3,
									    5,
									    3,
									    3,
									    0.15,
									    0.10,
									    0.98,
									    0.014,
									    0.035);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
							
			  TLorentzVector e1;
			  TLorentzVector e2;
			  TLorentzVector meson;
			  float mass = 0.;
			  for (unsigned int i=0; i<firstVector.size(); i++)
			    {
			      for (unsigned int j=0; j<secondVector.size() ; j++)
				{
									
				  if (firstVector[i] == secondVector[j])
				    continue;
				  e1.SetPtEtaPhiM(
						  ohPhotEt[firstVector[i]],
						  ohPhotEta[firstVector[i]],
						  ohPhotPhi[firstVector[i]],
						  0.);
				  e2.SetPtEtaPhiM(
						  ohPhotEt[secondVector[j]],
						  ohPhotEta[secondVector[j]],
						  ohPhotPhi[secondVector[j]],
						  0.);
				  meson = e1 + e2;
				  mass = meson.M();
									
				  if (mass>60)
				    triggerBit[it] = true;
									
				}
			    }
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoT_Photon_CaloIdL_IsoT) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    3,
									    5,
									    3,
									    3,
									    0.15,
									    0.10,
									    0.98,
									    0.014,
									    0.035);
		  if (firstVector.size()>=2)
		    {
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
      //warning: this is within if (triggerName.BeginsWith("OpenHLT_Photon26") )
      else if (triggerName.CompareTo(pathName_Photon_CaloIdL_IsoT_Photon) == 0)
	{
	  if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	    {
	      if (prescaleResponse(menu, cfg, rcounter, it))
		{
		  std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									    lowerEt,
									    0,
									    999,
									    999,
									    999,
									    999,
									    0.15,
									    0.10,
									    0.98);
		  if (firstVector.size()>=2)
		    {
		      //if(OpenHlt1PhotonPassed(upperEt,0,3,5,3,3,0.15,0.98,0.014,0.035)>=1) {      
		      std::vector<int> secondVector = VectorOpenHlt1PhotonPassed(
										 upperEt,
										 0,
										 3,
										 5,
										 3,
										 3,
										 0.15,
										 0.10,
										 0.98,
										 0.014,
										 0.035);
		      if (secondVector.size()>=1)
			{
			  triggerBit[it] = true;
			}
		    }
		}
	    }
	}
		
    }// end if(triggerName.BeginsWith("OpenHLT_Photon26") )
		
	
  // to be removed?
  // (added by Hartl when merging in new paths from V. Rekovic; most were obsolete but two using R9ID where new... including this)
  else if (triggerName.CompareTo("OpenHLT_Photon20_R9ID_Photon18_Isol_CaloId_L1R") == 0)
    {
      if (map_L1BitOfStandardHLTPath.find(triggerName)->second==1)
	{
	  if (prescaleResponse(menu, cfg, rcounter, it))
	    {
	      std::vector<int> firstVector = VectorOpenHlt1PhotonPassed(
									18,
									0,
									3,
									5,
									3,
									3,
									0.05,
									0.98,
									0.014,
									0.035);
	      if (firstVector.size()>=1)
		{
		  std::vector<int> secondVector =
		    VectorOpenHlt1PhotonPassedR9ID(
						   20,
						   0.8,
						   0,
						   999,
						   999,
						   999,
						   999,
						   0.15,
						   0.10,
						   0.98);
		  if (secondVector.size()>=1)
		    {
		      if ((firstVector.size() == 1)&&(secondVector.size()==1))
			{
			  if (firstVector.front() != secondVector.front())
			    triggerBit[it] = true;
			}
		      else
			triggerBit[it] = true;
		    }
		}
	    }
	}
    }
  else
    {
      if (nMissingTriggerWarnings < 100)
	cout << "Warning: the requested trigger " << triggerName
	     << " is not implemented in OHltTreeOpen. No rate will be calculated."
	     << endl;
      nMissingTriggerWarnings++;
    }
}

// functions


void OHltTree::PrintOhltVariables(int level, int type)
{
  cout << "Run " << Run <<", Event " << Event << endl;
  switch (type)
    {
    case muon:
			
      if (level == 3)
	{
				
	  cout << "Level 3: number of muons = " << NohMuL3 << endl;
				
	  for (int i=0; i<NohMuL3; i++)
	    {
	      cout << "ohMuL3Pt["<<i<<"] = " << ohMuL3Pt[i] << endl;
	      cout << "ohMuL3PtErr["<<i<<"] = " << ohMuL3PtErr[i] << endl;
	      cout << "ohMuL3Pt+Err["<<i<<"] = " << ohMuL3Pt[i]+2.2
		*ohMuL3PtErr[i]*ohMuL3Pt[i] << endl;
	      cout << "ohMuL3Phi["<<i<<"] = " << ohMuL3Phi[i] << endl;
	      cout << "ohMuL3Eta["<<i<<"] = " << ohMuL3Eta[i] << endl;
	      cout << "ohMuL3Chg["<<i<<"] = " << ohMuL3Chg[i] << endl;
	      cout << "ohMuL3Iso["<<i<<"] = " << ohMuL3Iso[i] << endl;
	      cout << "ohMuL3Dr["<<i<<"] = " << ohMuL3Dr[i] << endl;
	      cout << "ohMuL3Dz["<<i<<"] = " << ohMuL3Dz[i] << endl;
	      cout << "ohMuL3L2idx["<<i<<"] = " << ohMuL3L2idx[i] << endl;
	    }
	}
      else if (level == 2)
	{
	  cout << "Level 2: number of muons = " << NohMuL2 << endl;
	  for (int i=0; i<NohMuL2; i++)
	    {
	      cout << "ohMuL2Pt["<<i<<"] = " << ohMuL2Pt[i] << endl;
	      cout << "ohMuL2PtErr["<<i<<"] = " << ohMuL2PtErr[i] << endl;
	      cout << "ohMuL2Pt+Err["<<i<<"] = " << ohMuL2Pt[i]+3.9
		*ohMuL2PtErr[i]*ohMuL2Pt[i] << endl;
	      cout << "ohMuL2Phi["<<i<<"] = " << ohMuL2Phi[i] << endl;
	      cout << "ohMuL2Eta["<<i<<"] = " << ohMuL2Eta[i] << endl;
	      cout << "ohMuL2Chg["<<i<<"] = " << ohMuL2Chg[i] << endl;
	      cout << "ohMuL2Iso["<<i<<"] = " << ohMuL2Iso[i] << endl;
	      cout << "ohMuL2Dr["<<i<<"] = " << ohMuL2Dr[i] << endl;
	      cout << "ohMuL2Dz["<<i<<"] = " << ohMuL2Dz[i] << endl;
	    }
	}
      else if (level == 1)
	{
	  for (int i=0; i<NL1Mu; i++)
	    {
	      cout << "L1MuPt["<<i<<"] = " << L1MuPt[i] << endl;
	      cout << "L1MuEta["<<i<<"] = " << L1MuEta[i] << endl;
	      cout << "L1MuPhi["<<i<<"] = " << L1MuPhi[i] << endl;
	      cout << "L1MuIsol["<<i<<"] = " << L1MuIsol[i] << endl;
	      cout << "L1MuQal["<<i<<"] = " << L1MuQal[i] << endl;
	    }
	}
      else
	{
	  cout
	    << "PrintOhltVariables: Ohlt has Muon variables only for L1, 2, and 3. Must provide one."
	    << endl;
	}
      break;
			
    case electron:
      cout << "oh: number of electrons = " << NohEle << endl;
      for (int i=0; i<NohEle; i++)
	{
	  float ohElePt = ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1
									  *ohEleEta[i])));
	  cout << "ohEleEt["<<i<<"] = " << ohEleEt[i] << endl;
	  cout << "ohElePhi["<<i<<"] = " << ohElePhi[i] << endl;
	  cout << "ohEleEta["<<i<<"] = " << ohEleEta[i] << endl;
	  cout << "ohEleE["<<i<<"] = " << ohEleE[i] << endl;
	  cout << "ohEleP["<<i<<"] = " << ohEleP[i] << endl;
	  cout << "ohElePt["<<i<<"] =" << ohElePt << endl;
	  cout << "ohEleHiso["<<i<<"] = " << ohEleHiso[i] << endl;
	  cout << "ohEleTiso["<<i<<"] = " << ohEleTiso[i] << endl;
	  cout << "ohEleL1iso["<<i<<"] = " << ohEleL1iso[i] << endl;
	  cout << "ohEleHiso["<<i<<"]/ohEleEt["<<i<<"] = " << ohEleHiso[i]
	    /ohEleEt[i] << endl;
	  cout << "ohEleEiso["<<i<<"]/ohEleEt["<<i<<"] = " << ohEleEiso[i]
	    /ohEleEt[i] << endl;
	  cout << "ohEleTiso["<<i<<"]/ohEleEt["<<i<<"] = " << ohEleTiso[i]
	    /ohEleEt[i] << endl;
	  cout << "ohEleHforHoverE["<<i<<"] = " << ohEleHforHoverE[i] << endl;
	  cout << "ohEleHforHoverE["<<i<<"]/ohEleE["<<i<<"] = "
	       << ohEleHforHoverE[i]/ohEleE[i] << endl;
	  cout << "ohEleNewSC["<<i<<"] = " << ohEleNewSC[i] << endl;
	  cout << "ohElePixelSeeds["<<i<<"] = " << ohElePixelSeeds[i] << endl;
	  cout << "ohEleClusShap["<<i<<"] = " << ohEleClusShap[i] << endl;
	  cout << "ohEleR9["<<i<<"] = " << ohEleR9[i] << endl;
	  cout << "ohEleDeta["<<i<<"] = " << ohEleDeta[i] << endl;
	  cout << "ohEleDphi["<<i<<"] = " << ohEleDphi[i] << endl;
	}
			
      for (int i=0; i<NL1IsolEm; i++)
	{
	  cout << "L1IsolEmEt["<<i<<"] = " << L1IsolEmEt[i] << endl;
	  cout << "L1IsolEmE["<<i<<"] = " << L1IsolEmE[i] << endl;
	  cout << "L1IsolEmEta["<<i<<"] = " << L1IsolEmEta[i] << endl;
	  cout << "L1IsolEmPhi["<<i<<"] = " << L1IsolEmPhi[i] << endl;
	}
      for (int i=0; i<NL1NIsolEm; i++)
	{
	  cout << "L1NIsolEmEt["<<i<<"] = " << L1NIsolEmEt[i] << endl;
	  cout << "L1NIsolEmE["<<i<<"] = " << L1NIsolEmE[i] << endl;
	  cout << "L1NIsolEmEta["<<i<<"] = " << L1NIsolEmEta[i] << endl;
	  cout << "L1NIsolEmPhi["<<i<<"] = " << L1NIsolEmPhi[i] << endl;
	}
			
      break;
			
    case photon:
			
      cout << "oh: number of photons = " << NohPhot << endl;
      for (int i=0; i<NohPhot; i++)
	{
	  cout << "ohPhotEt["<<i<<"] = " << ohPhotEt[i] << endl;
	  cout << "ohPhotPhi["<<i<<"] = " << ohPhotPhi[i] << endl;
	  cout << "ohPhotEta["<<i<<"] = " << ohPhotEta[i] << endl;
	  cout << "ohPhotEiso["<<i<<"] = " << ohPhotEiso[i] << endl;
	  cout << "ohPhotHiso["<<i<<"] = " << ohPhotHiso[i] << endl;
	  cout << "ohPhotTiso["<<i<<"] = " << ohPhotTiso[i] << endl;
	  cout << "ohPhotL1iso["<<i<<"] = " << ohPhotL1iso[i] << endl;
	  cout << "ohPhotHiso["<<i<<"]/ohPhotEt["<<i<<"] = " << ohPhotHiso[i]
	    /ohPhotEt[i] << endl;
	  cout << "recoPhotE["<<i<<"] = " << recoPhotE[i] << endl;
	  cout << "recoPhotEt["<<i<<"] = " << recoPhotEt[i] << endl;
	  cout << "recoPhotPt["<<i<<"] = " << recoPhotPt[i] << endl;
	  cout << "recoPhotPhi["<<i<<"] = " << recoPhotPhi[i] << endl;
	  cout << "recoPhotEta["<<i<<"] = " << recoPhotEta[i] << endl;
				
	}
      break;
			
    case jet:
      cout << "oh: number of ohJetCal = " << NohJetCal << endl;
      for (int i=0; i<NohJetCal; i++)
	{
	  cout << "ohJetCalE["<<i<<"] = " << ohJetCalE[i] << endl;
	  cout << "ohJetCalPt["<<i<<"] = " << ohJetCalPt[i] << endl;
	  cout << "ohJetCalPhi["<<i<<"] = " << ohJetCalPhi[i] << endl;
	  cout << "ohJetCalEta["<<i<<"] = " << ohJetCalEta[i] << endl;
	}
      break;
			
    case tau:
      cout << "oh: number of taus = " << NohTau << endl;
      for (int i=0; i<NohTau; i++)
	{
	  cout<<"ohTauEt["<<i<<"] = " <<ohTauPt[i]<<endl;
	  cout<<"ohTauEiso["<<i<<"] = " <<ohTauEiso[i]<<endl;
	  cout<<"ohTauL25Tpt["<<i<<"] = " <<ohTauL25Tpt[i]<<endl;
	  cout<<"ohTauL25Tiso["<<i<<"] = " <<ohTauL25Tiso[i]<<endl;
	  cout<<"ohTauL3Tpt["<<i<<"] = " <<ohTauL3Tpt[i]<<endl;
	  cout<<"ohTauL3Tiso["<<i<<"] = " <<ohTauL3Tiso[i]<<endl;
	}
      break;
			
    default:
			
      cout << "PrintOhltVariables: You did not provide correct object type."
	   <<endl;
      break;
    }
}

// PFTau Leg
int OHltTree::OpenHltPFTauPassedNoMuon(
				       float Et,
				       float L25TrkPt,
				       float L3TrkIso,
				       float L3GammaIso)
{
  int rc = 0;
  // Loop over all oh ohpfTaus
  for (int i=0; i < NohpfTau; i++)
    {
      if (ohpfTauPt[i] >= Et && fabs(ohpfTauEta[i]) < 2.5
	  && ohpfTauLeadTrackPt[i] >= L25TrkPt && ohpfTauTrkIso[i] < L3TrkIso
	  && ohpfTauGammaIso[i] < L3GammaIso && OpenHltTauMuonMatching(
								       ohpfTauEta[i],
								       ohpfTauPhi[i]) == 0)
	rc++;
    }
  return rc;
}

int OHltTree::OpenHltPFTauPassedNoEle(
				      float Et,
				      float L25TrkPt,
				      int L3TrkIso,
				      int L3GammaIso)
{
	
  int rc = 0;
  // Loop over all oh ohpfTaus
  for (int i=0; i < NohpfTau; i++)
    {
      //        if (ohpfTauJetPt[i] >= Et) 
      if (ohpfTauPt[i] >= Et)
	if (fabs(ohpfTauEta[i])<2.5)
	  if (ohpfTauLeadTrackPt[i] >= L25TrkPt)
	    if (ohpfTauTrkIso[i] < L3TrkIso)
	      if (ohpfTauGammaIso[i] < L3GammaIso)
		if (OpenHltTauPFToCaloMatching(
					       ohpfTauEta[i],
					       ohpfTauPhi[i]) == 1)
		  if (OpenHltTauEleMatching(ohpfTauEta[i], ohpfTauPhi[i])
		      == 0)
		    rc++;
		
    }
	
  return rc;
}

int OHltTree::OpenHltTauMuonMatching(float eta, float phi)
{
  for (int j=0; j<NohMuL2; j++)
    {
      double deltaphi = fabs(phi-ohMuL2Phi[j]);
      if (deltaphi > 3.14159)
	deltaphi = (2.0 * 3.14159) - deltaphi;
      double deltaeta = fabs(eta-ohMuL2Eta[j]);
		
      if (sqrt(deltaeta*deltaeta + deltaphi*deltaphi) < 0.3)
	return 1;
    }
  return 0;
}

int OHltTree::OpenHltTauMuonMatching_wMuonID(float eta, float phi, double ptl1, double ptl2, double ptl3, double dr, int iso){
  // This example implements the new (CMSSW_2_X) flat muon pT cuts.
  // To emulate the old behavior, the cuts should be written
  // L2:        ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]
  // L3:        ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]
  int taumatch = 0;
  int rcL1 = 0; int rcL2 = 0; int rcL3 = 0; int rcL1L2L3 = 0;
  int NL1Mu = 8;
  int L1MinimalQuality = 4;
  int L1MaximalQuality = 7;
  int doL1L2matching = 0;
	
  for(int ic = 0; ic < 10; ic++)
    L3MuCandIDForOnia[ic] = -1;
	
  // Loop over all oh L3 muons and apply cuts
  for (int i=0;i<NohMuL3;i++) {
    int bestl1l2drmatchind = -1;
    double bestl1l2drmatch = 999.0;
		
    if( fabs(ohMuL3Eta[i]) < 2.5 ) { // L3 eta cut
      if(ohMuL3Pt[i] > ptl3) {  // L3 pT cut
	if(ohMuL3Dr[i] < dr) {  // L3 DR cut
	  if(ohMuL3Iso[i] >= iso) {  // L3 isolation
	    rcL3++;
						
	    // Begin L2 muons here.
	    // Get best L2<->L3 match, then
	    // begin applying cuts to L2
	    int j = ohMuL3L2idx[i];  // Get best L2<->L3 match
						
	    if ( (fabs(ohMuL2Eta[j])<2.5) ) {  // L2 eta cut
	      if( ohMuL2Pt[j] > ptl2 ) { // L2 pT cut
		if(ohMuL2Iso[j] >= iso) { // L2 isolation
		  rcL2++;
		  taumatch = j;
		  // Begin L1 muons here.
		  // Require there be an L1Extra muon Delta-R
		  // matched to the L2 candidate, and that it have
		  // good quality and pass nominal L1 pT cuts
		  for(int k = 0;k < NL1Mu;k++) {
		    if( (L1MuPt[k] < ptl1) ) // L1 pT cut
		      continue;
										
		    double deltaphi = fabs(ohMuL2Phi[j]-L1MuPhi[k]);
		    if(deltaphi > 3.14159)
		      deltaphi = (2.0 * 3.14159) - deltaphi;
										
		    double deltarl1l2 = sqrt((ohMuL2Eta[j]-L1MuEta[k])*(ohMuL2Eta[j]-L1MuEta[k]) +
					     (deltaphi*deltaphi));
		    if(deltarl1l2 < bestl1l2drmatch)
		      {
			bestl1l2drmatchind = k;
			bestl1l2drmatch = deltarl1l2;
		      }
		  } // End loop over L1Extra muons
									
		  if(doL1L2matching == 1)
		    {
		      // Cut on L1<->L2 matching and L1 quality
		      if((bestl1l2drmatch > 0.3) || (L1MuQal[bestl1l2drmatchind] < L1MinimalQuality) || (L1MuQal[bestl1l2drmatchind] > L1MaximalQuality))
			{
			  rcL1 = 0;
			  cout << "Failed L1-L2 match/quality" << endl;
			  cout << "L1-L2 delta-eta = " << L1MuEta[bestl1l2drmatchind] << ", " << ohMuL2Eta[j] << endl;
			  cout << "L1-L2 delta-pho = " << L1MuPhi[bestl1l2drmatchind] << ", " << ohMuL2Phi[j] << endl;
			  cout << "L1-L2 delta-R = " << bestl1l2drmatch << endl;
			}
		      else
			{
			  cout << "Passed L1-L2 match/quality" << endl;
			  L3MuCandIDForOnia[rcL1L2L3] = i;
			  rcL1++;
			  rcL1L2L3++;
			} // End L1 matching and quality cuts
		    }
		  else
		    {
		      L3MuCandIDForOnia[rcL1L2L3] = i;
		      rcL1L2L3++;
		    }
		} // End L2 isolation cut
	      } // End L2 eta cut
	    } // End L2 pT cut
	  } // End L3 isolation cut
	} // End L3 DR cut
      } // End L3 pT cut
    } // End L3 eta cut
  } // End loop over L3 muons
  double deltaphi = fabs(phi-ohMuL2Phi[taumatch]);
  if(deltaphi > 3.14159){ deltaphi = (2.0 * 3.14159) - deltaphi;}
  double deltaeta = fabs(eta-ohMuL2Eta[taumatch]);
  // if (deltaeta<0.3 && deltaphi<0.3){
  if(sqrt(deltaeta*deltaeta + deltaphi*deltaphi)<0.3){
    return 1;}
  return 0;
}


int OHltTree::OpenHltTauEleMatching(float eta, float phi)
{
  for (int j=0; j<NohEle; j++)
    {
      double deltaphi = fabs(phi-ohElePhi[j]);
      if (deltaphi > 3.14159)
	deltaphi = (2.0 * 3.14159) - deltaphi;
      double deltaeta = fabs(eta-ohEleEta[j]);
		
      if (deltaeta<0.3 && deltaphi<0.3)
	return 1;
    }
  return 0;
	
}

int OHltTree::OpenHltTauEleMatching_wEleID(float eta, float phi, float Et, int L1iso,
					   float Tisobarrel, float Tisoendcap,
					   float Tisoratiobarrel, float Tisoratioendcap,
					   float HisooverETbarrel, float HisooverETendcap,
					   float EisooverETbarrel, float EisooverETendcap,
					   float hoverebarrel, float hovereendcap,
					   float clusshapebarrel, float clusshapeendcap,
					   float r9barrel, float r9endcap,
					   float detabarrel, float detaendcap,
                                           float dphibarrel, float dphiendcap){
	
  {
    float barreleta = 1.479;
    float endcapeta = 2.65;
		
    // Loop over all oh electrons
    for (int i=0;i<NohEle;i++) {
      // ****************************************************
      // Bug fix
      // To be removed once the new ntuples are produced
      // ****************************************************
      float ohEleHoverE;
      float ohEleR9value;
      if(ohEleL1iso[i] == 1) {
	ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
	ohEleR9value = ohEleR9[i];
      }
      if(ohEleL1iso[i] == 0) {
	ohEleHoverE = ohEleR9[i]/ohEleE[i];
	ohEleR9value = ohEleHforHoverE[i];
      }
      // ****************************************************
      // ****************************************************
      int isbarrel = 0;
      int isendcap = 0;
      if(TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if(barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) < endcapeta)
	isendcap = 1;
			
			
      if ( ohEleEt[i] > Et) {
	if( TMath::Abs(ohEleEta[i]) < endcapeta ) {
	  if (ohEleNewSC[i]<=1) {
	    if (ohElePixelSeeds[i]>0) {
	      if ( ohEleL1iso[i] >= L1iso ) {  // L1iso is 0 or 1
		if( ohEleL1Dupl[i] == false) { // remove double-counted L1 SCs
		  if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETbarrel)) ||
		       (isendcap && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETendcap)) ) {
		    if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETbarrel)) ||
			 (isendcap && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETendcap)) ) {
		      if ( ((isbarrel) && (ohEleHoverE < hoverebarrel)) ||
			   ((isendcap) && (ohEleHoverE < hovereendcap))) {
			if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel && ohEleTiso[i] != -999.) || (Tisobarrel == 999.)))) ||
			     (isendcap && (((ohEleTiso[i] < Tisoendcap && ohEleTiso[i] != -999.) || (Tisoendcap == 999.))))) {
			  if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i] < Tisoratiobarrel)) ||
			      ((isendcap) && (ohEleTiso[i]/ohEleEt[i] < Tisoratioendcap))) {
			    if ( (isbarrel && ohEleClusShap[i] < clusshapebarrel) ||
				 (isendcap && ohEleClusShap[i] < clusshapeendcap) ) {
			      if ( (isbarrel && ohEleR9value < r9barrel) ||
				   (isendcap && ohEleR9value < r9endcap) ) {
				if ( (isbarrel && TMath::Abs(ohEleDeta[i]) < detabarrel) ||
				     (isendcap && TMath::Abs(ohEleDeta[i]) < detaendcap) ) {
				  if( (isbarrel && ohEleDphi[i] < dphibarrel) ||
				      (isendcap && ohEleDphi[i] < dphiendcap) ) {
				    double deltaphi = fabs(phi-ohElePhi[i]);
				    if(deltaphi > 3.14159){
				      deltaphi = (2.0 * 3.14159) - deltaphi;
				    }
				    double deltaeta = fabs(eta-ohEleEta[i]);
				    // if (deltaeta<0.3 && deltaphi<0.3){
				    if(sqrt(deltaeta*deltaeta + deltaphi*deltaphi)<0.3){
				      return 1;
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return 0;
}


int OHltTree::OpenHltTauPFToCaloMatching(float eta, float phi)
{
  for (int j=0; j<NohJetCal; j++)
    {
      if (ohJetCalPt[j]<8)
	continue;
      double deltaphi = fabs(phi-ohJetCalPhi[j]);
      if (deltaphi > 3.14159)
	deltaphi = (2.0 * 3.14159) - deltaphi;
      double deltaeta = fabs(eta-ohJetCalEta[j]);
		
      if (deltaeta<0.3 && deltaphi<0.3)
	return 1;
    }
  return 0;
	
}

int OHltTree::OpenHltL1L2TauMatching(
				     float eta,
				     float phi,
				     float tauThr,
				     float jetThr)
{
  for (int j=0; j<NL1Tau; j++)
    {
      double deltaphi = fabs(phi-L1TauPhi[j]);
      if (deltaphi > 3.14159)
	deltaphi = (2.0 * 3.14159) - deltaphi;
      double deltaeta = fabs(eta-L1TauEta[j]);
		
      if (deltaeta<0.3 && deltaphi<0.3 && L1TauEt[j]>tauThr)
	return 1;
    }
  for (int j=0; j<NL1CenJet; j++)
    {
      double deltaphi = fabs(phi-L1CenJetPhi[j]);
      if (deltaphi > 3.14159)
	deltaphi = (2.0 * 3.14159) - deltaphi;
      double deltaeta = fabs(eta-L1CenJetEta[j]);
		
      if (deltaeta<0.3 && deltaphi<0.3 && L1CenJetEt[j]>jetThr)
	return 1;
    }
  return 0;
}

int OHltTree::OpenHltTauPassed(
			       float Et,
			       float Eiso,
			       float L25Tpt,
			       int L25Tiso,
			       float L3Tpt,
			       int L3Tiso,
			       float L1TauThr,
			       float L1CenJetThr)
{
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohTau; i++)
    {
      if (ohTauPt[i] >= Et)
	{
	  if (ohTauEiso[i] <= Eiso)
	    if (ohTauL25Tpt[i] >= L25Tpt)
	      if (ohTauL25Tiso[i] >= L25Tiso)
		if (ohTauL3Tpt[i] >= L3Tpt)
		  if (ohTauL3Tiso[i] >= L3Tiso)
		    if (OpenHltL1L2TauMatching(
					       ohTauEta[i],
					       ohTauPhi[i],
					       L1TauThr,
					       L1CenJetThr) == 1)
		      rc++;
	}
    }
	
  return rc;
}

// L2 Ecal sliding cut isolation
int OHltTree::OpenHltTauL2SCPassed(
				   float Et,
				   float L25Tpt,
				   int L25Tiso,
				   float L3Tpt,
				   int L3Tiso,
				   float L1TauThr,
				   float L1CenJetThr)
{
  int rc = 0;
	
  // Loop over all oh electrons
  for (int i=0; i<NohTau; i++)
    {
      if (ohTauPt[i] >= Et)
	{
	  if (ohTauEiso[i] < (5 + 0.025*ohTauPt[i] + 0.0015*ohTauPt[i]
			      *ohTauPt[i])) // sliding cut
	    if (ohTauL25Tpt[i] >= L25Tpt)
	      if (ohTauL25Tiso[i] >= L25Tiso)
		if (ohTauL3Tpt[i] >= L3Tpt)
		  if (ohTauL3Tiso[i] >= L3Tiso)
		    if (OpenHltL1L2TauMatching(
					       ohTauEta[i],
					       ohTauPhi[i],
					       L1TauThr,
					       L1CenJetThr) == 1)
		      rc++;
	}
    }
	
  return rc;
}


//NVV
int OHltTree::OpenHltIsoPFTauPassed(float Et, float eta, float LTpT, float L1TauThr, float L1CenJetThr){
  int count = 0;

  // Isolation thresholds
  float maxTrkIso		= 1.0;
  float maxGammaIso	= 1.5;

  // Loop over all ohpfTaus
  for (int i = 0; i < NohpfTau; i++){
    if (ohpfTauPt[i] >= Et){ // min ET
      if (fabs(ohpfTauEta[i]) < eta){ // eta constraint
	if( (ohpfTauTrkIso[i] < maxTrkIso) && (ohpfTauGammaIso[i] < maxGammaIso) ){ // Isolation (track and gamma)
	  if (ohpfTauLeadTrackPt[i] >= LTpT){ // min leading track pT
	    if (OpenHltL1L2TauMatching( ohpfTauEta[i], ohpfTauPhi[i], L1TauThr, L1CenJetThr) == 1){ // Matching to L1 object
	      count++;
	    }
	  }
	}
      }
    }
  }

  return count;
}

int OHltTree::OpenHltTightConeIsoPFTauPassed(float Et, float eta, float LTpT, float L1TauThr, float L1CenJetThr){
  int count = 0;

  // Isolation thresholds
  float maxTrkIso		= 1.0;
  float maxGammaIso	= 1.5;

  // Loop over all ohpfTaus
  for (int i = 0; i < NohpfTau; i++){
    if (ohpfTauPt[i] >= Et){ // min ET
      if (fabs(ohpfTauEta[i]) < eta){ // eta constraint
	if( (ohpfTauTightConeTrkIso[i] < maxTrkIso) && (ohpfTauTightConeGammaIso[i] < maxGammaIso) ){ // Isolation (track and gamma)
	  if (ohpfTauLeadTrackPt[i] >= LTpT){ // min leading track pT
	    if (OpenHltL1L2TauMatching( ohpfTauEta[i], ohpfTauPhi[i], L1TauThr, L1CenJetThr) == 1){ // Matching to L1 object
	      count++;
	    }
	  }
	}
      }
    }
  }

  return count;
}



int OHltTree::OpenHltTauL2SCMETPassed(
				      float Et,
				      float L25Tpt,
				      int L25Tiso,
				      float L3Tpt,
				      int L3Tiso,
				      float met,
				      float L1TauThr,
				      float L1CenJetThr)
{
  int rc = 0;
	
  // Loop over all oh electrons
  for (int i=0; i<NohTau; i++)
    {
      if (ohTauPt[i]>= Et)
	{
	  if (ohTauEiso[i]
	      < (5 + 0.025*ohTauPt[i] + 0.0015*ohTauPt[i]*ohTauPt[i])) // sliding cut
	    if (ohTauL25Tpt[i]>= L25Tpt)
	      if (ohTauL25Tiso[i]>= L25Tiso)
		if (ohTauL3Tpt[i]>= L3Tpt)
		  if (ohTauL3Tiso[i]>= L3Tiso)
		    if (OpenHltL1L2TauMatching(
					       ohTauEta[i],
					       ohTauPhi[i],
					       L1TauThr,
					       L1CenJetThr) == 1)
		      if (recoMetCal> met)
			rc++;
	}
    }
  return rc;
}

int OHltTree::OpenHlt2Tau1LegL3IsoPassed(
					 float Et,
					 float L25Tpt,
					 int L25Tiso,
					 float L3Tpt,
					 float L1TauThr,
					 float L1CenJetThr)
{
  int rc = 0;
  int l3iso = 0;
	
  // Loop over all oh taus
  for (int i=0; i<NohTau; i++)
    {
      if (ohTauPt[i] >= Et)
	{
	  if (ohTauEiso[i] < (5 + 0.025*ohTauPt[i] + 0.0015*ohTauPt[i]
			      *ohTauPt[i])) // sliding cut
	    if (ohTauL25Tpt[i] >= L25Tpt)
	      if (ohTauL25Tiso[i] >= L25Tiso)
		if (ohTauL3Tpt[i] >= L3Tpt)
		  if (OpenHltL1L2TauMatching(
					     ohTauEta[i],
					     ohTauPhi[i],
					     L1TauThr,
					     L1CenJetThr) == 1)
		    {
		      rc++;
		      if (ohTauL3Tiso[i] >= 1)
			l3iso++;
		    }
	}
    }
	
  if (rc>=2)
    return l3iso;
  return 0;
}

int OHltTree::OpenHlt1Ele1PFTauPassed(
				      float Et, float Eta, int L1iso,
				      float Tisobarrel, float Tisoendcap,
				      float Tisoratiobarrel, float Tisoratioendcap,
				      float HisooverETbarrel, float HisooverETendcap,
				      float EisooverETbarrel, float EisooverETendcap,
				      float hoverebarrel, float hovereendcap,
				      float clusshapebarrel, float clusshapeendcap,
				      float r9barrel, float r9endcap,
				      float detabarrel, float detaendcap,
				      float dphibarrel, float dphiendcap,
				      float TauEt, float TauEta, float L25TrkPt,
				      float L3TrkIso, float L3GammaIso, float PFMHTCut)
{
  int rc=0;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      bool isBarrel, isEndCap, doTaus;
      isBarrel = isEndCap = doTaus = false;
		
      if (fabs(ohEleEta[i]) < 1.479)
	isBarrel = true;
      if (fabs(ohEleEta[i]) > 1.479 && fabs(ohEleEta[i]) < 2.65)
	isEndCap = true;
		
      if (ohEleEt[i] > Et && 
	  fabs(ohEleEta[i]) < Eta && 
	  (isBarrel || isEndCap) && 
	  ohEleNewSC[i] == 1 &&
	  ohElePixelSeeds[i] > 0 &&
	  ohEleL1iso[i] >= L1iso && // L1iso is 0 or 1
	  !ohEleL1Dupl[i]) { // remove double-counted L1 SCs
	if (isBarrel) {
	  if (ohEleHiso[i]/ohEleEt[i] < HisooverETbarrel &&
	      ohEleEiso[i]/ohEleEt[i] < EisooverETbarrel && 
	      ohEleHoverE < hoverebarrel &&
	      ((ohEleTiso[i] < Tisobarrel && ohEleTiso[i] != -999.) || Tisobarrel == 999.) && 
	      ohEleTiso[i]/ohEleEt[i] < Tisoratiobarrel && 
	      ohEleClusShap[i] < clusshapebarrel && 
	      ohEleR9[i] < r9barrel && 
	      fabs(ohEleDeta[i]) < detabarrel && 
	      ohEleDphi[i] < dphibarrel)
	    doTaus = true;
	}
	if (isEndCap) {
	  if (ohEleHiso[i]/ohEleEt[i] < HisooverETendcap && 
	      ohEleEiso[i]/ohEleEt[i] < EisooverETendcap && 
	      ohEleHoverE < hovereendcap && 
	      ((ohEleTiso[i] < Tisoendcap && ohEleTiso[i] != -999.) || (Tisoendcap == 999.)) && 
	      ohEleTiso[i]/ohEleEt[i] < Tisoratioendcap && 
	      ohEleClusShap[i] < clusshapeendcap &&
	      ohEleR9[i] < r9endcap && fabs(ohEleDeta[i]) < detabarrel && 
	      ohEleDphi[i] < detaendcap)
	    doTaus = true;
	}
      }
      if (doTaus)
	{
	  for (int j=0; j < NohpfTau; j++)
	    {
	      if (ohpfTauPt[j] >= TauEt && 
		  fabs(ohpfTauEta[j]) < TauEta && 
		  ohpfTauLeadTrackPt[j] >= L25TrkPt && 
		  ohpfTauTrkIso[j] < L3TrkIso && 
		  ohpfTauGammaIso[j] < L3GammaIso && 
		  pfMHT >= PFMHTCut)
		{
		  float dphi = fabs(ohElePhi[i] - ohpfTauPhi[j]);
		  float deta = fabs(ohEleEta[i] - ohpfTauEta[j]);
		  if (dphi > 3.14159)
		    dphi = (2.0 * 3.14159) - dphi;
		  if (sqrt(dphi*dphi + deta*deta) > 0.3)
		    ++rc; 
		}
	    }
	}
    }
  return rc;
}

int OHltTree::OpenHlt1PhotonSamHarperPassed(
					    float Et,
					    int L1iso,
					    float Tisobarrel,
					    float Tisoendcap,
					    float Tisoratiobarrel,
					    float Tisoratioendcap,
					    float HisooverETbarrel,
					    float HisooverETendcap,
					    float EisooverETbarrel,
					    float EisooverETendcap,
					    float hoverebarrel,
					    float hovereendcap,
					    float clusshapebarrel,
					    float clusshapeendcap,
					    float r9barrel,
					    float r9endcap,
					    float detabarrel,
					    float detaendcap,
					    float dphibarrel,
					    float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      float quadraticEcalIsol = ohPhotEiso[i] - (0.012 * ohPhotEt[i]);
      float quadraticHcalIsol = ohPhotHiso[i] - (0.005 * ohPhotEt[i]);
      float quadraticTrackIsol = ohPhotTiso[i] - (0.002 * ohPhotEt[i]);
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < endcapeta)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  if (ohPhotL1Dupl[i] == false)
		    { // remove double-counted L1 SCs 
		      if ( (isbarrel && (quadraticHcalIsol < HisooverETbarrel))
			   || (isendcap && (quadraticHcalIsol < HisooverETendcap)))
			{
			  if ( (isbarrel && (quadraticEcalIsol < EisooverETbarrel))
			       || (isendcap && (quadraticEcalIsol
						< EisooverETendcap)))
			    {
			      if ( ((isbarrel) && (ohPhotHoverE < hoverebarrel))
				   || ((isendcap) && (ohPhotHoverE < hovereendcap)))
				{
				  if (((isbarrel) && (quadraticTrackIsol < Tisobarrel))
				      || ((isendcap) && (quadraticTrackIsol
							 < Tisoendcap)))
				    {
				      if ( (isbarrel && ohPhotClusShap[i]
					    < clusshapebarrel) || (isendcap
								   && ohPhotClusShap[i] < clusshapeendcap))
					{
					  if ( (isbarrel && ohPhotR9[i] < r9barrel)
					       || (isendcap && ohPhotR9[i] < r9endcap))
					    {
					      rc++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

//Lucie
int OHltTree::OpenHlt1PhotonPassed(
				   float Et,
				   std::map< TString, float> r9Id,
				   std::map< TString, float> caloId,
				   std::map< TString, float> photonIso
				   )
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
  int L1iso = 0;
	
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      float quadraticEcalIsol = ohPhotEiso[i] - (0.012 * ohPhotEt[i]);
      float quadraticHcalIsol = ohPhotHiso[i] - (0.005 * ohPhotEt[i]);
      float quadraticTrackIsol = ohPhotTiso[i] - (0.002 * ohPhotEt[i]);
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < endcapeta)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  if (ohPhotL1Dupl[i] == false)
		    { // remove double-counted L1 SCs 
		      if ( (isbarrel && (quadraticHcalIsol < photonIso["HisoBR"]))
			   || (isendcap && (quadraticHcalIsol < photonIso["HisoEC"] )))
			{
			  if ( (isbarrel && (quadraticEcalIsol <  photonIso["Eiso"]))
			       || (isendcap && (quadraticEcalIsol
						< photonIso["Eiso"])))
			    {
			      if ( ((isbarrel) && (ohPhotHoverE < caloId["hoverebarrel"]))
				   || ((isendcap) && (ohPhotHoverE < caloId["hovereendcap"])))
				{
				  if (((isbarrel) && (quadraticTrackIsol < photonIso["Tiso"]))
				      || ((isendcap) && (quadraticTrackIsol
							 < photonIso["Tiso"] )))
				    {
				      if ( (isbarrel && ohPhotClusShap[i]
					    < caloId["clusshapebarrel"]) || (isendcap
								   && ohPhotClusShap[i] < caloId["clusshapeendcap"]))
					{
					  if ( (isbarrel && ohPhotR9[i] < r9Id["HoverEEB"])
					       || (isendcap && ohPhotR9[i] < r9Id["HoverEEC"]))
					    {
					      rc++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}


//Lucie / Arnaud 
vector<int> OHltTree::VecOpenHlt1PhotonPassed(
				   float Et,
				   std::map< TString, float> r9Id,
				   std::map< TString, float> caloId,
				   std::map< TString, float> photonIso
				   )
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
  int L1iso = 0;
	
  vector<int> rc;

  // Loop over all oh electrons
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      float quadraticEcalIsol = ohPhotEiso[i] - (0.012 * ohPhotEt[i]);
      float quadraticHcalIsol = ohPhotHiso[i] - (0.005 * ohPhotEt[i]);
      float quadraticTrackIsol = ohPhotTiso[i] - (0.002 * ohPhotEt[i]);
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < endcapeta)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  if (ohPhotL1Dupl[i] == false)
		    { // remove double-counted L1 SCs 
		      if ( (isbarrel && (quadraticHcalIsol < photonIso["HisoBR"]))
			   || (isendcap && (quadraticHcalIsol < photonIso["HisoEC"] )))
			{
			  if ( (isbarrel && (quadraticEcalIsol <  photonIso["Eiso"]))
			       || (isendcap && (quadraticEcalIsol
						< photonIso["Eiso"])))
			    {
			      if ( ((isbarrel) && (ohPhotHoverE < caloId["hoverebarrel"]))
				   || ((isendcap) && (ohPhotHoverE < caloId["hovereendcap"])))
				{
				  if (((isbarrel) && (quadraticTrackIsol < photonIso["Tiso"]))
				      || ((isendcap) && (quadraticTrackIsol
							 < photonIso["Tiso"] )))
				    {
				      if ( (isbarrel && ohPhotClusShap[i]
					    < caloId["clusshapebarrel"]) || (isendcap
								   && ohPhotClusShap[i] < caloId["clusshapeendcap"]))
					{
					  if ( (isbarrel && ohPhotR9[i] < r9Id["HoverEEB"])
					       || (isendcap && ohPhotR9[i] < r9Id["HoverEEC"]))
					    {
					      rc.push_back(i);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}



//Arnaud ohEcalActiv
vector<int>  OHltTree::VecOpenHlt1EcalActivPassed(
				   float Et,
				   std::map< TString, float> r9Id,
				   std::map< TString, float> caloId,
				   std::map< TString, float> photonIso
				   )
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
  int L1iso = 0;
	
  vector<int> rc;

  // Loop over all oh electrons
  for (int i=0; i<NohEcalActiv; i++)
    {
      float ohEcalActivE = ohEcalActivEt[i] / (sin(2*atan(exp(-1.0*ohEcalActivEta[i]))));
      float ohEcalActivHoverE = ohEcalActivHforHoverE[i]/ohEcalActivE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEcalActivEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEcalActivEta[i]) && TMath::Abs(ohEcalActivEta[i])
	  < endcapeta)
	isendcap = 1;
		
      float quadraticEcalIsol = ohEcalActivEiso[i] - (0.012 * ohEcalActivEt[i]);
      float quadraticHcalIsol = ohEcalActivHiso[i] - (0.005 * ohEcalActivEt[i]);
      float quadraticTrackIsol = ohEcalActivTiso[i] - (0.002 * ohEcalActivEt[i]);
		
      if (ohEcalActivEt[i] > Et)
	{
	  if (TMath::Abs(ohEcalActivEta[i]) < endcapeta)
	    {
	      if (ohEcalActivL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  //if (ohEcalActivL1Dupl[i] == false ){
		     // remove double-counted L1 SCs 
		      if ( (isbarrel && (quadraticHcalIsol < photonIso["HisoBR"]))

			   || (isendcap && (quadraticHcalIsol < photonIso["HisoEC"] )))
			{
			  if ( (isbarrel && (quadraticEcalIsol <  photonIso["Eiso"]))
			       || (isendcap && (quadraticEcalIsol
						< photonIso["Eiso"])))
			    {
			      if ( ((isbarrel) && (ohEcalActivHoverE < caloId["hoverebarrel"]))
				   || ((isendcap) && (ohEcalActivHoverE < caloId["hovereendcap"])))
				{
				  if (((isbarrel) && (quadraticTrackIsol < photonIso["Tiso"]))
				      || ((isendcap) && (quadraticTrackIsol
							 < photonIso["Tiso"] )))
				    {
				      if ( (isbarrel && ohEcalActivClusShap[i]
					    < caloId["clusshapebarrel"]) || (isendcap
								   && ohEcalActivClusShap[i] < caloId["clusshapeendcap"]))
					{
					  if ( (isbarrel && ohEcalActivR9[i] < r9Id["HoverEEB"])
					       || (isendcap && ohEcalActivR9[i] < r9Id["HoverEEC"]))
					    {
					      rc.push_back(i);
					    }
					}
				    }
				}
			    }
			}
		      //}
		}
	    }
	}
    }
	
  return rc;
}

vector<int> OHltTree::VectorOpenHlt1PhotonSamHarperPassed(
							  float Et,
							  int L1iso,
							  float Tisobarrel,
							  float Tisoendcap,
							  float Tisoratiobarrel,
							  float Tisoratioendcap,
							  float HisooverETbarrel,
							  float HisooverETendcap,
							  float EisooverETbarrel,
							  float EisooverETendcap,
							  float hoverebarrel,
							  float hovereendcap,
							  float clusshapebarrel,
							  float clusshapeendcap,
							  float r9barrel,
							  float r9endcap,
							  float detabarrel,
							  float detaendcap,
							  float dphibarrel,
							  float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  vector<int> rc;
  // Loop over all oh electrons
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      float quadraticEcalIsol = ohPhotEiso[i] - (0.012 * ohPhotEt[i]);
      float quadraticHcalIsol = ohPhotHiso[i] - (0.005 * ohPhotEt[i]);
      float quadraticTrackIsol = ohPhotTiso[i] - (0.002 * ohPhotEt[i]);
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < endcapeta)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  if (ohPhotL1Dupl[i] == false)
		    { // remove double-counted L1 SCs 
		      if ( (isbarrel && (quadraticHcalIsol < HisooverETbarrel))
			   || (isendcap && (quadraticHcalIsol < HisooverETendcap)))
			{
			  if ( (isbarrel && (quadraticEcalIsol < EisooverETbarrel))
			       || (isendcap && (quadraticEcalIsol
						< EisooverETendcap)))
			    {
			      if ( ((isbarrel) && (ohPhotHoverE < hoverebarrel))
				   || ((isendcap) && (ohPhotHoverE < hovereendcap)))
				{
				  if (((isbarrel) && (quadraticTrackIsol < Tisobarrel))
				      || ((isendcap) && (quadraticTrackIsol
							 < Tisoendcap)))
				    {
				      if ( (isbarrel && ohPhotClusShap[i]
					    < clusshapebarrel) || (isendcap
								   && ohPhotClusShap[i] < clusshapeendcap))
					{
					  if ( (isbarrel && ohPhotR9[i] < r9barrel)
					       || (isendcap && ohPhotR9[i] < r9endcap))
					    {
					      rc.push_back(i);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHltPhoCuts(
			     float e_t,
			     float he_eb,
			     float he_ee,
			     float shape_eb,
			     float shape_ee,
			     float cal_iso, 
			     float trk_iso)
{
	
  // Default cleaning cuts for all photon paths
  float r9barrel = 999.;//0.98;
  float r9endcap = 999.;//1.0;
  float EtaGap = 1.479;
  float EtaMax = 2.65;
	
  int rc = 0;
	
  // Loop over all oh photons
  for (int i=0; i<NohPhot; i++)
    {
      if (ohPhotL1Dupl[i])
	continue;
		
      float eta = TMath::Abs(ohPhotEta[i]);
      if (eta > EtaMax)
	continue;
		
      if (ohPhotEt[i] < e_t)
	continue;
		
      bool isBarrel = (eta < EtaGap);
      bool isEndcap = (eta >= EtaGap);
		
      bool passSpikeCleaning = (isBarrel && ohPhotR9[i] < r9barrel)
	|| (isEndcap && ohPhotR9[i] < r9endcap);
      if (!passSpikeCleaning)
	continue;
		
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
		
      bool passHoverE = (isBarrel && ohPhotHoverE < he_eb) || (isEndcap
							       && ohPhotHoverE < he_ee);
      if (!passHoverE)
	continue;
		
      bool passShape = (isBarrel && ohPhotClusShap[i] < shape_eb) || (isEndcap
								      && ohPhotClusShap[i] < shape_ee);
      if (!passShape)
	continue;
		
      if (ohPhotEiso[i]/ohPhotEt[i] > cal_iso)
	continue;
      if (ohPhotHiso[i]/ohPhotEt[i] > cal_iso)
	continue;
      if (ohPhotTiso[i]/ohPhotEt[i] > trk_iso)
	continue;
		
      rc++;
    }// for
	
  return rc;
}

int OHltTree::OpenHlt1PhotonPassedRA3(
				      float Et,
				      int L1iso,
				      float HisooverETbarrel,
				      float HisooverETendcap,
				      float EisooverETbarrel,
				      float EisooverETendcap,
				      float hoverebarrel,
				      float hovereendcap,
				      float r9barrel,
				      float r9endcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  int rc = 0;
  // Loop over all oh photons
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < endcapeta)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{ // L1iso is 0 or 1 
		  if (ohPhotL1Dupl[i] == false)
		    { // remove double-counted L1 SCs 
		      if ( (isbarrel && ((ohPhotHiso[i]/ohPhotEt[i])
					 < HisooverETbarrel)) || (isendcap && ((ohPhotHiso[i]
										/ohPhotEt[i]) < HisooverETendcap)))
			{
			  if ( (isbarrel && ((ohPhotEiso[i]/ohPhotEt[i])
					     < EisooverETbarrel)) || (isendcap && ((ohPhotEiso[i]
										    /ohPhotEt[i]) < EisooverETendcap)))
			    {
			      if ( ((isbarrel) && (ohPhotHoverE < hoverebarrel))
				   || ((isendcap) && (ohPhotHoverE < hovereendcap)))
				{
				  if ( (isbarrel && ohPhotR9[i] < r9barrel)
				       || (isendcap && ohPhotR9[i] < r9endcap))
				    {
				      rc++;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}

bool OHltTree::OpenHLT_EleX_R9cut(const float& Et, const float& r9barrel)
{
  return OpenHlt1ElectronSamHarperPassed(Et, 0, // ET, L1isolation 
					 999.,
					 999., // Track iso barrel, Track iso endcap 
					 999.,
					 999., // Track/pT iso barrel, Track/pT iso endcap 
					 999.,
					 999., // H/ET iso barrel, H/ET iso endcap 
					 999.,
					 999., // E/ET iso barrel, E/ET iso endcap 
					 0.15,
					 0.15, // H/E barrel, H/E endcap 
					 999.,
					 999., // cluster shape barrel, cluster shape endcap 
					 r9barrel,
					 1.0, // R9 barrel, R9 endcap 
					 999.,
					 999., // Deta barrel, Deta endcap 
					 999.,
					 999. // Dphi barrel, Dphi endcap
					 ) >= 1;
}



int OHltTree::OpenHlt1ElectronSamHarperPassed(
					      float Et,
					      int L1iso,
					      float Tisobarrel,
					      float Tisoendcap,
					      float Tisoratiobarrel,
					      float Tisoratioendcap,
					      float HisooverETbarrel,
					      float HisooverETendcap,
					      float EisooverETbarrel,
					      float EisooverETendcap,
					      float hoverebarrel,
					      float hovereendcap,
					      float clusshapebarrel,
					      float clusshapeendcap,
					      float r9barrel,
					      float r9endcap,
					      float detabarrel,
					      float detaendcap,
					      float dphibarrel,
					      float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]<=1)
		{
		  if (ohElePixelSeeds[i]>0)
		    {
		      if (ohEleL1iso[i] >= L1iso)
			{ // L1iso is 0 or 1 
			  if (ohEleL1Dupl[i] == false)
			    { // remove double-counted L1 SCs 
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
						     < EisooverETbarrel)) || (isendcap
									      && ((ohEleEiso[i]/ohEleEt[i])
										  < EisooverETendcap)))
				    {
				      if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
					   || ((isendcap) && (ohEleHoverE
							      < hovereendcap)))
					{
					  if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel
							       && ohEleTiso[i] != -999.) || (Tisobarrel
											     == 999.)))) || (isendcap
													     && (((ohEleTiso[i] < Tisoendcap
														   && ohEleTiso[i] != -999.)
														  || (Tisoendcap == 999.)))))
					    {
					      if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i]
								  < Tisoratiobarrel)) || ((isendcap)
											  && (ohEleTiso[i]/ohEleEt[i]
											      < Tisoratioendcap)))
						{
						  if ( (isbarrel && ohEleClusShap[i]
							< clusshapebarrel) || (isendcap
									       && ohEleClusShap[i]
									       < clusshapeendcap))
						    {
						      if ( (isbarrel && ohEleR9[i]
							    < r9barrel) || (isendcap
									    && ohEleR9[i] < r9endcap))
							{
							  if ( (isbarrel
								&& TMath::Abs(ohEleDeta[i])
								< detabarrel)
							       || (isendcap
								   && TMath::Abs(ohEleDeta[i])
								   < detaendcap))
							    {
							      if ( (isbarrel && ohEleDphi[i]
								    < dphibarrel)
								   || (isendcap
								       && ohEleDphi[i]
								       < dphiendcap))
								{
								  rc++;
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

//Lucie
int OHltTree::OpenHlt1ElectronPassed(float Et,
				     std::map< TString, float> caloId,
				     std::map< TString, float> caloIso,
				     std::map< TString, float> trkId,
				     std::map< TString, float> trkIso
				     )
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
  float Tisobarrel = 999.;
  float Tisoendcap = 999.;
  int L1iso = 0;	

  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
	
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]<=1)
		{
		  if (ohElePixelSeeds[i]>0)
		    {
		      if (ohEleL1iso[i] >= L1iso)
			{ // L1iso is 0 or 1 
			  if (ohEleL1Dupl[i] == false)
			    { // remove double-counted L1 SCs 
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < caloIso["HisooverETbarrel"]))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < caloIso["HisooverETendcap"])))
				{
				  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
						     < caloIso["EisooverETbarrel"])) || (isendcap
									      && ((ohEleEiso[i]/ohEleEt[i])
										  < caloIso["EisooverETendcap"])))
				    {
				      if ( ((isbarrel) && (ohEleHoverE < caloId["hoverebarrel"]))
					   || ((isendcap) && (ohEleHoverE
							      < caloId["hovereendcap"])))
					{
					  if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel
							       && ohEleTiso[i] != -999.) || (Tisobarrel
											     == 999.)))) || (isendcap
													     && (((ohEleTiso[i] < Tisoendcap
														   && ohEleTiso[i] != -999.)
														  || (Tisoendcap == 999.)))))
					    {
					      if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i]
								  < trkIso["Tisoratiobarrel"])) || ((isendcap)
											  && (ohEleTiso[i]/ohEleEt[i]
											      < trkIso["Tisoratioendcap"])))
						{
						  if ( (isbarrel && ohEleClusShap[i]
							< caloId["clusshapebarrel"]) || (isendcap
									       && ohEleClusShap[i]
									       < caloId["clusshapeendcap"]))
						    {
						      if ( isbarrel || isendcap)
							{
							  if ( (isbarrel
								&& TMath::Abs(ohEleDeta[i])
								< trkId["detabarrel"])
							       || (isendcap
								   && TMath::Abs(ohEleDeta[i])
								   < trkId["detaendcap"]))
							    {
							      if ( (isbarrel && ohEleDphi[i]
								    < trkId["dphibarrel"])
								   || (isendcap
								       && ohEleDphi[i]
								       < trkId["dphiendcap"]))
								{
								  rc++;
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}


int OHltTree::OpenHlt2ElectronsSamHarperPassed(
					       float Et,
					       int L1iso,
					       float Tisobarrel,
					       float Tisoendcap,
					       float Tisoratiobarrel,
					       float Tisoratioendcap,
					       float HisooverETbarrel,
					       float HisooverETendcap,
					       float EisooverETbarrel,
					       float EisooverETendcap,
					       float hoverebarrel,
					       float hovereendcap,
					       float clusshapebarrel,
					       float clusshapeendcap,
					       float r9barrel,
					       float r9endcap,
					       float detabarrel,
					       float detaendcap,
					       float dphibarrel,
					       float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  int rc = 0;
  int rcsconly = 0;
	
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]==1)
		{
		  if (ohEleL1iso[i] >= L1iso)
		    { // L1iso is 0 or 1 
		      if (ohEleL1Dupl[i] == false)
			{ // remove double-counted L1 SCs 
			  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
					     < EisooverETbarrel)) || (isendcap && ((ohEleEiso[i]
										    /ohEleEt[i]) < EisooverETendcap)))
			    {
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
				       || ((isendcap) && (ohEleHoverE < hovereendcap)))
				    {
				      if ( (isbarrel && ohEleClusShap[i]
					    < clusshapebarrel) || (isendcap
								   && ohEleClusShap[i] < clusshapeendcap))
					{
					  if ( (isbarrel && ohEleR9[i] < r9barrel)
					       || (isendcap && ohEleR9[i] < r9endcap))
					    {
					      if (ohElePixelSeeds[i]>0)
						{
						  rcsconly++;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  if (rcsconly >= 2)
    {
      for (int i=0; i<NohEle; i++)
	{
	  float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
	  int isbarrel = 0;
	  int isendcap = 0;
	  if (TMath::Abs(ohEleEta[i]) < barreleta)
	    isbarrel = 1;
	  if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	      < endcapeta)
	    isendcap = 1;
			
	  if (ohEleEt[i] > Et)
	    {
	      if (TMath::Abs(ohEleEta[i]) < endcapeta)
		{
		  if (ohEleNewSC[i]<=1)
		    {
		      if (ohElePixelSeeds[i]>0)
			{
			  if (ohEleL1iso[i] >= L1iso)
			    { // L1iso is 0 or 1 
			      if (ohEleL1Dupl[i] == false)
				{ // remove double-counted L1 SCs 
				  if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						     < HisooverETbarrel)) || (isendcap
									      && ((ohEleHiso[i]/ohEleEt[i])
										  < HisooverETendcap)))
				    {
				      if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
							 < EisooverETbarrel)) || (isendcap
										  && ((ohEleEiso[i]/ohEleEt[i])
										      < EisooverETendcap)))
					{
					  if ( ((isbarrel) && (ohEleHoverE
							       < hoverebarrel)) || ((isendcap)
										    && (ohEleHoverE < hovereendcap)))
					    {
					      if ( (isbarrel
						    && (((ohEleTiso[i] < Tisobarrel
							  && ohEleTiso[i] != -999.)
							 || (Tisobarrel == 999.))))
						   || (isendcap && (((ohEleTiso[i]
								      < Tisoendcap && ohEleTiso[i]
								      != -999.) || (Tisoendcap
										    == 999.)))))
						{
						  if (((isbarrel) && (ohEleTiso[i]
								      /ohEleEt[i] < Tisoratiobarrel))
						      || ((isendcap) && (ohEleTiso[i]
									 /ohEleEt[i]
									 < Tisoratioendcap)))
						    {
						      if ( (isbarrel && ohEleClusShap[i]
							    < clusshapebarrel) || (isendcap
										   && ohEleClusShap[i]
										   < clusshapeendcap))
							{
							  if ( (isbarrel && ohEleR9[i]
								< r9barrel) || (isendcap
										&& ohEleR9[i] < r9endcap))
							    {
							      if ( (isbarrel
								    && TMath::Abs(ohEleDeta[i])
								    < detabarrel)
								   || (isendcap
								       && TMath::Abs(ohEleDeta[i])
								       < detaendcap))
								{
								  if ( (isbarrel
									&& ohEleDphi[i]
									< dphibarrel)
								       || (isendcap
									   && ohEleDphi[i]
									   < dphiendcap))
								    {
								      rc++;
								    }
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt3ElectronsSamHarperPassed(
					       float Et,
					       int L1iso,
					       float Tisobarrel,
					       float Tisoendcap,
					       float Tisoratiobarrel,
					       float Tisoratioendcap,
					       float HisooverETbarrel,
					       float HisooverETendcap,
					       float EisooverETbarrel,
					       float EisooverETendcap,
					       float hoverebarrel,
					       float hovereendcap,
					       float clusshapebarrel,
					       float clusshapeendcap,
					       float r9barrel,
					       float r9endcap,
					       float detabarrel,
					       float detaendcap,
					       float dphibarrel,
					       float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  int rc = 0;
  int rcsconly = 0;
	
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      // bug fix
      float ohEleHoverE;
      float ohEleR9v;
      //if(ohEleL1iso[i] == 1){
      ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      ohEleR9v = ohEleR9[i];
      //}
		
      //    float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]==1)
		{
		  if (ohEleL1iso[i] >= L1iso)
		    { // L1iso is 0 or 1 
		      if (ohEleL1Dupl[i] == false)
			{ // remove double-counted L1 SCs 
			  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
					     < EisooverETbarrel)) || (isendcap && ((ohEleEiso[i]
										    /ohEleEt[i]) < EisooverETendcap)))
			    {
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
				       || ((isendcap) && (ohEleHoverE < hovereendcap)))
				    {
				      if ( (isbarrel && ohEleClusShap[i]
					    < clusshapebarrel) || (isendcap
								   && ohEleClusShap[i] < clusshapeendcap))
					{
					  //                      if ( (isbarrel && ohEleR9[i] < r9barrel) ||
					  //                         (isendcap && ohEleR9[i] < r9endcap) ) {
					  if ( (isbarrel && ohEleR9v < r9barrel)
					       || (isendcap && ohEleR9v < r9endcap))
					    {
					      if (ohElePixelSeeds[i]>0)
						{
						  rcsconly++;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  if (rcsconly >= 3)
    {
      for (int i=0; i<NohEle; i++)
	{
	  // bug fix
	  float ohEleHoverE;
	  float ohEleR9v;
	  //if(ohEleL1iso[i] == 1){    
	  ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
	  ohEleR9v = ohEleR9[i];
			
	  int isbarrel = 0;
	  int isendcap = 0;
	  if (TMath::Abs(ohEleEta[i]) < barreleta)
	    isbarrel = 1;
	  if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	      < endcapeta)
	    isendcap = 1;
	  if (ohEleEt[i] > Et)
	    {
	      if (TMath::Abs(ohEleEta[i]) < endcapeta)
		{
		  if (ohEleNewSC[i]<=1)
		    {
		      if (ohElePixelSeeds[i]>0)
			{
			  if (ohEleL1iso[i] >= L1iso)
			    { // L1iso is 0 or 1 
			      if (ohEleL1Dupl[i] == false)
				{ // remove double-counted L1 SCs 
				  if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						     < HisooverETbarrel)) || (isendcap
									      && ((ohEleHiso[i]/ohEleEt[i])
										  < HisooverETendcap)))
				    {
				      if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
							 < EisooverETbarrel)) || (isendcap
										  && ((ohEleEiso[i]/ohEleEt[i])
										      < EisooverETendcap)))
					{
					  if ( ((isbarrel) && (ohEleHoverE
							       < hoverebarrel)) || ((isendcap)
										    && (ohEleHoverE < hovereendcap)))
					    {
					      if ( (isbarrel
						    && (((ohEleTiso[i] < Tisobarrel
							  && ohEleTiso[i] != -999.)
							 || (Tisobarrel == 999.))))
						   || (isendcap && (((ohEleTiso[i]
								      < Tisoendcap && ohEleTiso[i]
								      != -999.) || (Tisoendcap
										    == 999.)))))
						{
						  if (((isbarrel) && (ohEleTiso[i]
								      /ohEleEt[i] < Tisoratiobarrel))
						      || ((isendcap) && (ohEleTiso[i]
									 /ohEleEt[i]
									 < Tisoratioendcap)))
						    {
						      if ( (isbarrel && ohEleClusShap[i]
							    < clusshapebarrel) || (isendcap
										   && ohEleClusShap[i]
										   < clusshapeendcap))
							{
							  //                              if ( (isbarrel && ohEleR9[i] < r9barrel) ||
							  //                                 (isendcap && ohEleR9[i] < r9endcap) ) {
							  if ( (isbarrel && ohEleR9v
								< r9barrel) || (isendcap
										&& ohEleR9v < r9endcap))
							    {
							      if ( (isbarrel
								    && TMath::Abs(ohEleDeta[i])
								    < detabarrel)
								   || (isendcap
								       && TMath::Abs(ohEleDeta[i])
								       < detaendcap))
								{
								  if ( (isbarrel
									&& ohEleDphi[i]
									< dphibarrel)
								       || (isendcap
									   && ohEleDphi[i]
									   < dphiendcap))
								    {
								      rc++;
								    }
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHltGetElectronsSamHarperPassed(
						 int *Passed,
						 float Et,
						 int L1iso,
						 float Tisobarrel,
						 float Tisoendcap,
						 float Tisoratiobarrel,
						 float Tisoratioendcap,
						 float HisooverETbarrel,
						 float HisooverETendcap,
						 float EisooverETbarrel,
						 float EisooverETendcap,
						 float hoverebarrel,
						 float hovereendcap,
						 float clusshapebarrel,
						 float clusshapeendcap,
						 float r9barrel,
						 float r9endcap,
						 float detabarrel,
						 float detaendcap,
						 float dphibarrel,
						 float dphiendcap)
{
  int NPassed = 0;
	
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  // First check if only one electron is going to pass the cuts using the one electron code
  // we use the one electron code to look to see if 0 or 1 electrons pass this condition
  int rc = 0;
  int tmpindex=-999;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      //    float ohEleHoverE = ohEleHiso[i]/ohEleEt[i];
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]<=1)
		{
		  if (ohElePixelSeeds[i]>0)
		    {
		      if (ohEleL1iso[i] >= L1iso)
			{ // L1iso is 0 or 1 
			  if (ohEleL1Dupl[i] == false)
			    { // remove double-counted L1 SCs 
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
						     < EisooverETbarrel)) || (isendcap
									      && ((ohEleEiso[i]/ohEleEt[i])
										  < EisooverETendcap)))
				    {
				      if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
					   || ((isendcap) && (ohEleHoverE
							      < hovereendcap)))
					{
					  if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel
							       && ohEleTiso[i] != -999.) || (Tisobarrel
											     == 999.)))) || (isendcap
													     && (((ohEleTiso[i] < Tisoendcap
														   && ohEleTiso[i] != -999.)
														  || (Tisoendcap == 999.)))))
					    {
					      if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i]
								  < Tisoratiobarrel)) || ((isendcap)
											  && (ohEleTiso[i]/ohEleEt[i]
											      < Tisoratioendcap)))
						{
						  if ( (isbarrel && ohEleClusShap[i]
							< clusshapebarrel) || (isendcap
									       && ohEleClusShap[i]
									       < clusshapeendcap))
						    {
						      if ( (isbarrel && ohEleR9[i]
							    < r9barrel) || (isendcap
									    && ohEleR9[i] < r9endcap))
							{
							  if ( (isbarrel
								&& TMath::Abs(ohEleDeta[i])
								< detabarrel)
							       || (isendcap
								   && TMath::Abs(ohEleDeta[i])
								   < detaendcap))
							    {
							      if ( (isbarrel && ohEleDphi[i]
								    < dphibarrel)
								   || (isendcap
								       && ohEleDphi[i]
								       < dphiendcap))
								{
								  rc++;
								  tmpindex=i; // temporarily store the index of this event
								  if (rc >1)
								    break; // If we have 2 electrons passing, we need the 2 ele code
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  //  cout << "rc = " << rc << endl;;
  if (rc == 0)
    {
      return 0;
    }
  if (rc == 1)
    {
      Passed[NPassed++]=tmpindex; // if only 1 ele matched we can use this result without looping on the 2 ele code
      return NPassed;
    }
  // otherwise, we use the 2 ele code:
	
  int rcsconly=0;
  std::vector<int> csPassedEle;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      //float ohEleHoverE = ohEleHiso[i]/ohEleEt[i];
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]==1)
		{
		  if (ohEleL1iso[i] >= L1iso)
		    { // L1iso is 0 or 1 
		      if (ohEleL1Dupl[i] == false)
			{ // remove double-counted L1 SCs 
			  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
					     < EisooverETbarrel)) || (isendcap && ((ohEleEiso[i]
										    /ohEleEt[i]) < EisooverETendcap)))
			    {
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
				       || ((isendcap) && (ohEleHoverE < hovereendcap)))
				    {
				      if ( (isbarrel && ohEleClusShap[i]
					    < clusshapebarrel) || (isendcap
								   && ohEleClusShap[i] < clusshapeendcap))
					{
					  if ( (isbarrel && ohEleR9[i] < r9barrel)
					       || (isendcap && ohEleR9[i] < r9endcap))
					    {
					      if (ohElePixelSeeds[i]>0)
						{
						  rcsconly++;
						  csPassedEle.push_back(i);
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  //  cout << "rcsconly = " << rcsconly << endl;
  if (rcsconly == 0)
    { // This really shouldn't happen, but included for safety
      return NPassed;
    }
  if (rcsconly == 1)
    { // ok, we only had 1 cs, but 2 eles were assigned to it
      Passed[NPassed++] = tmpindex;
      return NPassed;
    }
	
  if (rcsconly >= 2)
    {
      for (int i=0; i<NohEle; i++)
	{
	  float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
	  int isbarrel = 0;
	  int isendcap = 0;
	  if (TMath::Abs(ohEleEta[i]) < barreleta)
	    isbarrel = 1;
	  if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	      < endcapeta)
	    isendcap = 1;
			
	  if (ohEleEt[i] > Et)
	    {
	      if (TMath::Abs(ohEleEta[i]) < endcapeta)
		{
		  if (ohEleNewSC[i]<=1)
		    {
		      if (ohElePixelSeeds[i]>0)
			{
			  if (ohEleL1iso[i] >= L1iso)
			    { // L1iso is 0 or 1 
			      if (ohEleL1Dupl[i] == false)
				{ // remove double-counted L1 SCs 
				  if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						     < HisooverETbarrel)) || (isendcap
									      && ((ohEleHiso[i]/ohEleEt[i])
										  < HisooverETendcap)))
				    {
				      if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
							 < EisooverETbarrel)) || (isendcap
										  && ((ohEleEiso[i]/ohEleEt[i])
										      < EisooverETendcap)))
					{
					  if ( ((isbarrel) && (ohEleHoverE
							       < hoverebarrel)) || ((isendcap)
										    && (ohEleHoverE < hovereendcap)))
					    {
					      if ( (isbarrel
						    && (((ohEleTiso[i] < Tisobarrel
							  && ohEleTiso[i] != -999.)
							 || (Tisobarrel == 999.))))
						   || (isendcap && (((ohEleTiso[i]
								      < Tisoendcap && ohEleTiso[i]
								      != -999.) || (Tisoendcap
										    == 999.)))))
						{
						  if (((isbarrel) && (ohEleTiso[i]
								      /ohEleEt[i] < Tisoratiobarrel))
						      || ((isendcap) && (ohEleTiso[i]
									 /ohEleEt[i]
									 < Tisoratioendcap)))
						    {
						      if ( (isbarrel && ohEleClusShap[i]
							    < clusshapebarrel) || (isendcap
										   && ohEleClusShap[i]
										   < clusshapeendcap))
							{
							  if ( (isbarrel && ohEleR9[i]
								< r9barrel) || (isendcap
										&& ohEleR9[i] < r9endcap))
							    {
							      if ( (isbarrel
								    && TMath::Abs(ohEleDeta[i])
								    < detabarrel)
								   || (isendcap
								       && TMath::Abs(ohEleDeta[i])
								       < detaendcap))
								{
								  if ( (isbarrel
									&& ohEleDphi[i]
									< dphibarrel)
								       || (isendcap
									   && ohEleDphi[i]
									   < dphiendcap))
								    {
								      for (unsigned int j=0; j
									     <csPassedEle.size(); j++)
									{
									  if (i
									      == csPassedEle.at(j))
									    { // check if the electron is in the cs matching list
									      Passed[NPassed++]
										= i;
									      rc++; // ok, don't really need this, but keeping for debugging
									      break;
									    }
									}
								    }
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return NPassed;
}

int OHltTree::OpenHlt2ElectronsAsymSamHarperPassed(
						   float Et1,
						   int L1iso1,
						   float Tisobarrel1,
						   float Tisoendcap1,
						   float Tisoratiobarrel1,
						   float Tisoratioendcap1,
						   float HisooverETbarrel1,
						   float HisooverETendcap1,
						   float EisooverETbarrel1,
						   float EisooverETendcap1,
						   float hoverebarrel1,
						   float hovereendcap1,
						   float clusshapebarrel1,
						   float clusshapeendcap1,
						   float r9barrel1,
						   float r9endcap1,
						   float detabarrel1,
						   float detaendcap1,
						   float dphibarrel1,
						   float dphiendcap1,
						   float Et2,
						   int L1iso2,
						   float Tisobarrel2,
						   float Tisoendcap2,
						   float Tisoratiobarrel2,
						   float Tisoratioendcap2,
						   float HisooverETbarrel2,
						   float HisooverETendcap2,
						   float EisooverETbarrel2,
						   float EisooverETendcap2,
						   float hoverebarrel2,
						   float hovereendcap2,
						   float clusshapebarrel2,
						   float clusshapeendcap2,
						   float r9barrel2,
						   float r9endcap2,
						   float detabarrel2,
						   float detaendcap2,
						   float dphibarrel2,
						   float dphiendcap2)
{
  //  cout << "AB" << endl;
  int FirstEle[8000], SecondEle[8000];
  // cout << "BA" << endl;  
  int NFirst = OpenHltGetElectronsSamHarperPassed(
						  FirstEle,
						  Et1,
						  L1iso1,
						  Tisobarrel1,
						  Tisoendcap1,
						  Tisoratiobarrel1,
						  Tisoratioendcap1,
						  HisooverETbarrel1,
						  HisooverETendcap1,
						  EisooverETbarrel1,
						  EisooverETendcap1,
						  hoverebarrel1,
						  hovereendcap1,
						  clusshapebarrel1,
						  clusshapeendcap1,
						  r9barrel1,
						  r9endcap1,
						  detabarrel1,
						  detaendcap1,
						  dphibarrel1,
						  dphiendcap1);
  int NSecond = OpenHltGetElectronsSamHarperPassed(
						   SecondEle,
						   Et2,
						   L1iso2,
						   Tisobarrel2,
						   Tisoendcap2,
						   Tisoratiobarrel2,
						   Tisoratioendcap2,
						   HisooverETbarrel2,
						   HisooverETendcap2,
						   EisooverETbarrel2,
						   EisooverETendcap2,
						   hoverebarrel2,
						   hovereendcap2,
						   clusshapebarrel2,
						   clusshapeendcap2,
						   r9barrel2,
						   r9endcap2,
						   detabarrel2,
						   detaendcap2,
						   dphibarrel2,
						   dphiendcap2);
  // std::cout << "ABBA " << NFirst << "  " << NSecond << endl;
  if (NFirst == 0 || NSecond == 0)
    return 0; // if no eles passed one condition, fail
  if (NFirst == 1 && NSecond == 1 && FirstEle[0] == SecondEle[0]) 
    return 0; //only 1 electron passed both conditions
  return 1; // in any other case, at least 1 unique electron passed each condition, so pass the event
	
}

int OHltTree::OpenHlt1ElectronPassed(float Et, int L1iso, float Tiso, float Hiso)
{
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < 2.65)
	    //	if ( ohEleHiso[i] < Hiso || ohEleHiso[i]/ohEleEt[i] < 0.05)
	    if (ohEleHiso[i]/ohEleEt[i] < 0.15)
	      if (ohEleNewSC[i]==1)
		if (ohElePixelSeeds[i]>0)
		  if ( (ohEleTiso[i] < Tiso && ohEleTiso[i] != -999.)
		       || (Tiso == 9999.))
		    if (ohEleL1iso[i] >= L1iso) // L1iso is 0 or 1
		      if (ohEleL1Dupl[i] == false) // remove double-counted L1 SCs  
			rc++;
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1PhotonPassed(
				   float Et,
				   int L1iso,
				   float Tiso,
				   float Eiso,
				   float HisoBR,
				   float HisoEC)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  // Default cleaning cuts for all photon paths
  float r9barrel = 0.98;
  float r9endcap = 1.0;
  float hoverebarrel = 0.15;
  float hovereendcap = 0.15;
  int rc = 0;
	
  // Loop over all oh photons 
  for (int i=0; i<NohPhot; i++)
    {
      float ohPhotE = ohPhotEt[i] / (sin(2*atan(exp(-1.0*ohPhotEta[i]))));
      float ohPhotHoverE = ohPhotHforHoverE[i]/ohPhotE;
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohPhotEt[i] > Et)
	{
	  if ( ((isbarrel) && (ohPhotHoverE < hoverebarrel)) || ((isendcap)
								 && (ohPhotHoverE < hovereendcap)))
	    {
	      if ( (isbarrel && ohPhotR9[i] < r9barrel) || (isendcap
							    && ohPhotR9[i] < r9endcap))
		{
		  if (ohPhotL1Dupl[i] == false) // remove double-counted L1 SCs    
		    rc++;
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt2ElectronMassWinPassed(
					    float Et,
					    int L1iso,
					    float Hiso,
					    float massLow,
					    float massHigh)
{
  TLorentzVector e1;
  TLorentzVector e2;
  TLorentzVector meson;
	
  int rc = 0;
	
  for (int i=0; i<NohEle; i++)
    {
      for (int j=0; j<NohEle && j != i; j++)
	{
	  if (ohEleEt[i] > Et && ohEleEt[j] > Et)
	    {
	      if (TMath::Abs(ohEleEta[i]) < 2.65 && TMath::Abs(ohEleEta[j])
		  < 2.65)
		{
		  if ( ((ohEleHiso[i] < Hiso) || (ohEleHiso[i]/ohEleEt[i] < 0.2))
		       && ((ohEleHiso[j] < Hiso) || (ohEleHiso[j]/ohEleEt[j]
						     < 0.2)))
		    {
		      if (ohEleNewSC[i]==1 && ohEleNewSC[j]==1)
			{
			  if (ohElePixelSeeds[i]>0 && ohElePixelSeeds[j]>0)
			    {
			      if (ohEleL1iso[i] >= L1iso && ohEleL1iso[j] >= L1iso)
				{ // L1iso is 0 or 1  
				  if (ohEleL1Dupl[i] == false && ohEleL1Dupl[j]
				      == false)
				    { // remove double-counted L1 SCs    
				      e1.SetPtEtaPhiM(
						      ohEleEt[i],
						      ohEleEta[i],
						      ohElePhi[i],
						      0.0);
				      e2.SetPtEtaPhiM(
						      ohEleEt[j],
						      ohEleEta[j],
						      ohElePhi[j],
						      0.0);
				      meson = e1 + e2;
										
				      float mesonmass = meson.M();
				      if (mesonmass > massLow && mesonmass < massHigh)
					rc++;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  if (rc>0)
    return 1;
  else
    return 0;
}

int OHltTree::OpenHlt2Electron1LegIdPassed(
					   float Et,
					   int L1iso,
					   float Tiso,
					   float Hiso)
{
  int rc = 0;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < 2.65)
	    {
	      if (ohEleHiso[i] < Hiso || ohEleHiso[i]/ohEleEt[i] < 0.05)
		{
		  if (ohEleNewSC[i]==1)
		    {
		      if (ohElePixelSeeds[i]>0)
			{
			  if ( (ohEleTiso[i] < Tiso && ohEleTiso[i] != -999.)
			       || (Tiso == 9999.))
			    {
			      if (ohEleL1iso[i] >= L1iso)
				{ // L1iso is 0 or 1
				  if (ohEleL1Dupl[i] == false)
				    { // remove double-counted L1 SCs
				      // Loop over all oh electrons
				      for (int i=0; i<NohEle; i++)
					{
					  if (ohEleEt[i] > Et)
					    {
					      if (TMath::Abs(ohEleEta[i]) < 2.65)
						if (ohEleHiso[i] < Hiso || ohEleHiso[i]
						    /ohEleEt[i] < 0.05)
						  if (ohEleNewSC[i]==1)
						    if (ohElePixelSeeds[i]>0)
						      if ( (ohEleTiso[i] < Tiso
							    && ohEleTiso[i] != -999.)
							   || (Tiso == 9999.))
							if (ohEleL1iso[i] >= L1iso) // L1iso is 0 or 1
							  if ( (TMath::Abs(ohEleEta[i])
								< 1.479
								&& ohEleClusShap[i]
								< 0.015)
							       || (1.479
								   < TMath::Abs(ohEleEta[i])
								   && TMath::Abs(ohEleEta[i])
								   < 2.65
								   && ohEleClusShap[i]
								   < 0.04))
							    if ( (ohEleDeta[i]
								  < 0.008)
								 && (ohEleDphi[i]
								     < 0.1))
							      if (ohEleL1Dupl[i]
								  == false) // remove double-counted L1 SCs
								rc++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1ElectronHTPassed(
				       float Et,
				       float HT,
				       float jetThreshold,
				       int L1iso,
				       float Tiso,
				       float Hiso,
				       float dr)
{
  vector<int> PassedElectrons;
  int NPassedEle=0;
  for (int i=0; i<NohEle; i++)
    {
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < 2.65)
	    //if ( ohEleHiso[i] < Hiso || ohEleHiso[i]/ohEleEt[i] < 0.05)
	    if (ohEleHiso[i]/ohEleEt[i] < 0.15)
	      if (ohEleNewSC[i]==1)
		if (ohElePixelSeeds[i]>0)
		  if ( (ohEleTiso[i] < Tiso && ohEleTiso[i] != -999.)
		       || (Tiso == 9999.))
		    if (ohEleL1iso[i] >= L1iso) // L1iso is 0 or 1
		      if (ohEleL1Dupl[i] == false)
			{ // remove double-counted L1 SCs  
			  PassedElectrons.push_back(i);
			  NPassedEle++;
			}
			
	}
    }
  if (NPassedEle==0)
    return 0;
	
  float sumHT=0;
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i] < jetThreshold)
	continue;
      bool jetPass=true;
      for (unsigned int iEle = 0; iEle<PassedElectrons.size(); iEle++)
	{
	  float dphi = ohElePhi[PassedElectrons.at(iEle)] - ohJetCalPhi[i];
	  float deta = ohEleEta[PassedElectrons.at(iEle)] - ohJetCalEta[i];
	  if (dphi*dphi+deta*deta<dr*dr) // require electron not in any jet
	    jetPass=false;
	}
      if (jetPass)
	sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
    }
  if (sumHT>HT)
    return 1;
  return 0;
}

int OHltTree::OpenHlt1ElectronEleIDHTPassed(
					    float Et,
					    float HT,
					    float jetThreshold,
					    int L1iso,
					    float Tiso,
					    float Hiso,
					    float dr)
{
  vector<int> PassedElectrons;
  int NPassedEle=0;
  for (int i=0; i<NohEle; i++)
    {
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < 2.65)
	    //if ( ohEleHiso[i] < Hiso || ohEleHiso[i]/ohEleEt[i] < 0.05)
	    if (ohEleHiso[i]/ohEleEt[i] < 0.15)
	      if (ohEleNewSC[i]==1)
		if (ohElePixelSeeds[i]>0)
		  if ( (ohEleTiso[i] < Tiso && ohEleTiso[i] != -999.)
		       || (Tiso == 9999.))
		    if (ohEleL1iso[i] >= L1iso) // L1iso is 0 or 1
		      if ( (TMath::Abs(ohEleEta[i]) < 1.479
			    && ohEleClusShap[i] < 0.015) || (1.479
							     < TMath::Abs(ohEleEta[i])
							     && TMath::Abs(ohEleEta[i]) < 2.65
							     && ohEleClusShap[i] < 0.04))
			if ( (ohEleDeta[i] < 0.008) && (ohEleDphi[i]
							< 0.1))
			  if (ohEleL1Dupl[i] == false)
			    { // remove double-counted L1 SCs  
			      PassedElectrons.push_back(i);
			      NPassedEle++;
			    }
			
	}
    }
  if (NPassedEle==0)
    return 0;
	
  float sumHT=0;
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i] < jetThreshold)
	continue;
      bool jetPass=true;
      for (unsigned int iEle = 0; iEle<PassedElectrons.size(); iEle++)
	{
	  float dphi = ohElePhi[PassedElectrons.at(iEle)] - ohJetCalPhi[i];
	  float deta = ohEleEta[PassedElectrons.at(iEle)] - ohJetCalEta[i];
	  if (dphi*dphi+deta*deta<dr*dr) // require electron not in any jet
	    jetPass=false;
	}
      if (jetPass)
	sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
    }
  if (sumHT>HT)
    return 1;
  return 0;
}

int OHltTree::OpenHlt1BJetPassedEleRemoval(
      float jetEt,
      float jetEta,
      float drcut,
      float discL25,
      float discL3,
      float Et,
      int L1iso,
      float Tisobarrel,
      float Tisoendcap,
      float Tisoratiobarrel,
      float Tisoratioendcap,
      float HisooverETbarrel,
      float HisooverETendcap,
      float EisooverETbarrel,
      float EisooverETendcap,
      float hoverebarrel,
      float hovereendcap,
      float clusshapebarrel,
      float clusshapeendcap,
      float r9barrel,
      float r9endcap,
      float detabarrel,
      float detaendcap,
      float dphibarrel,
      float dphiendcap)
{

  int rc = 0;
  int njets = 0;
  
   //Loop over corrected oh b-jets
   for (int j = 0; j < NohBJetL2Corrected; j++)
   {//loop over jets
     
     bool isOverlapping = false;
     
     // ****************************************************
     // Exclude jets which are matched to electrons
     // ****************************************************
     float barreleta = 1.479;
     float endcapeta = 2.65;
     
     // Loop over all oh electrons
     for (int i=0; i<NohEle; i++)
       {//loop over electrons

	 int isbarrel = 0;
	 int isendcap = 0;
	 if (TMath::Abs(ohEleEta[i]) < barreleta)
	   isbarrel = 1;
	 if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	     < endcapeta)
	   isendcap = 1;
	 
	 if (ohEleEt[i] > Et)
	   {
	     if (TMath::Abs(ohEleEta[i]) < endcapeta)
               {
		 if (ohEleNewSC[i]<=1)
		   {
                     if (ohElePixelSeeds[i]>0)
		       {
			 if (ohEleL1iso[i] >= L1iso)
			   { // L1iso is 0 or 1 
			     if (ohEleL1Dupl[i] == false)
			       { // remove double-counted L1 SCs 
				 if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETbarrel)) || 
				      (isendcap && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETendcap)))
				   {
				     if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETbarrel)) || 
					  (isendcap && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETendcap)))
				       {
					 if ( ((isbarrel) && (ohEleHforHoverE[i]/ohEleE[i] < hoverebarrel)) || 
					      ((isendcap) && (ohEleHforHoverE[i]/ohEleE[i] < hovereendcap)))
					   {
					     if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel && ohEleTiso[i] != -999.) || 
								 (Tisobarrel == 999.))))
						  || (isendcap && (((ohEleTiso[i] < Tisoendcap && ohEleTiso[i] != -999.) || 
								    (Tisoendcap == 999.)))))
					       {
						 if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i] < Tisoratiobarrel)) || 
						     ((isendcap) && (ohEleTiso[i]/ohEleEt[i] < Tisoratioendcap)))
						   {
						     if ( (isbarrel && ohEleClusShap[i] < clusshapebarrel)
							  || (isendcap && ohEleClusShap[i] < clusshapeendcap))
						       {
							 if ( (isbarrel && ohEleR9[i] < r9barrel) || 
							      (isendcap && ohEleR9[i] < r9endcap))
							   {
							     if ( (isbarrel && TMath::Abs(ohEleDeta[i]) < detabarrel) || 
								  (isendcap && TMath::Abs(ohEleDeta[i]) < detaendcap))
							       {
								 if ( (isbarrel && ohEleDphi[i] < dphibarrel) || 
								      (isendcap && ohEleDphi[i] < dphiendcap))
								   {
								     
								     double deltaphi = fabs(ohBJetL2CorrectedPhi[j] - ohElePhi[i]);
								     if (deltaphi > 3.14159)
								       deltaphi = (2.0 * 3.14159) - deltaphi;
								     
								     double deltaRJetEle = sqrt((ohBJetL2CorrectedEta[j]-ohEleEta[i])
												*(ohBJetL2CorrectedEta[j]-ohEleEta[i])
												+ (deltaphi*deltaphi));
								     
								     if (deltaRJetEle < drcut)
								       {
									 isOverlapping = true;
									 break;
								       }
								   }
							       }
							   }
						       }
						   }
					       }
					   }
				       }
				   }
			       }
			   }
		       }
		   }
               }
	   }
       }//loop over electrons
     
     if (!isOverlapping)
       {//overlap
	 
	 if (ohBJetL2CorrectedEt[j] > jetEt && fabs(ohBJetL2CorrectedEta[j]) < jetEta) {// ET and eta cuts
	   
	   if (ohBJetIPL25Tag[j] >= discL25)
	     { // Level 2.5 b tag  
	       if (ohBJetIPL3Tag[j] >= discL3)
		 { // Level 3 b tag  
		   njets++;
		 }
	     }
	 }//ET and eta cuts
       }//overlap  
   }//loop over jets
   
   if (njets >= 1) rc = true;

   return rc;
}

bool OHltTree::OpenHltNCorJetPassedEleRemoval(
      int N,
      float jetPt,
      float jetEta,
      float drcut,
      float Et,
      int L1iso,
      float Tisobarrel,
      float Tisoendcap,
      float Tisoratiobarrel,
      float Tisoratioendcap,
      float HisooverETbarrel,
      float HisooverETendcap,
      float EisooverETbarrel,
      float EisooverETendcap,
      float hoverebarrel,
      float hovereendcap,
      float clusshapebarrel,
      float clusshapeendcap,
      float r9barrel,
      float r9endcap,
      float detabarrel,
      float detaendcap,
      float dphibarrel,
      float dphiendcap)
{

   float barreleta = 1.479;
   float endcapeta = 2.65;

   // Loop over all oh electrons
   for (int i=0; i<NohEle; i++)
      {
        int isbarrel = 0;
        int isendcap = 0;
        if (TMath::Abs(ohEleEta[i]) < barreleta)
            isbarrel = 1;
        if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) < endcapeta)
               isendcap = 1;

        if (ohEleEt[i] > Et)
            {
               if (TMath::Abs(ohEleEta[i]) < endcapeta)
               {
                  if (ohEleNewSC[i]<=1)
                  {
                     if (ohElePixelSeeds[i]>0)
                     {
                        if (ohEleL1iso[i] >= L1iso)
                        { // L1iso is 0 or 1 
                           if (ohEleL1Dupl[i] == false)
                           { // remove double-counted L1 SCs 
                              if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETbarrel)) ||
                                   (isendcap && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETendcap)))
                              {
                                 if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETbarrel)) ||
                                    (isendcap  && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETendcap)))
                                 {
                                    if ( ((isbarrel) && (ohEleHforHoverE[i]/ohEleE[i] < hoverebarrel)) ||
                                    	 ((isendcap) && (ohEleHforHoverE[i]/ohEleE[i] < hovereendcap)))
                                    {
                                       if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel && ohEleTiso[i] != -999.) ||
                                                           (Tisobarrel == 999.))))
                                             || (isendcap && (((ohEleTiso[i] < Tisoendcap && ohEleTiso[i] != -999.) ||
                                                               (Tisoendcap  == 999.)))))
                                       {
                                          if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i] < Tisoratiobarrel)) ||
                                              ((isendcap) && (ohEleTiso[i]/ohEleEt[i] < Tisoratioendcap)))
                                          {
                                             if ( (isbarrel && ohEleClusShap[i]  < clusshapebarrel)
                                                   || (isendcap && ohEleClusShap[i] < clusshapeendcap))
                                             {
                                                if ( (isbarrel && ohEleR9[i] < r9barrel) ||
                                                     (isendcap && ohEleR9[i] < r9endcap))
                                                {
                                                   if ( (isbarrel && TMath::Abs(ohEleDeta[i]) < detabarrel) ||
                                                        (isendcap && TMath::Abs(ohEleDeta[i]) < detaendcap))
                                                   {
                                                      if ( (isbarrel && ohEleDphi[i] < dphibarrel) ||
                                                           (isendcap && ohEleDphi[i] < dphiendcap))
                                                      {
							
							int nGoodJets = 0; 
							//Loop over all oh corrected jets
							for (int j = 0; j < NohJetCorCal; j++)
							{
							  if (ohJetCorCalPt[j] > jetPt && fabs(ohJetCorCalEta[j]) < jetEta && OpenJetID(j))
							    { // PT, eta and JetID cuts
							       double deltaphi = fabs(ohJetCorCalPhi[j] -ohElePhi[i]);
							       if (deltaphi > 3.14159) deltaphi = (2.0 * 3.14159) - deltaphi;
							       double deltaRJetEle = sqrt((ohJetCorCalEta[j]-ohEleEta[i])
                                                                                   *(ohJetCorCalEta[j]-ohEleEta[i])
                                                                                   + (deltaphi *deltaphi));

							       if (!(deltaRJetEle < drcut)) nGoodJets++;
							    }
							}  
							if (nGoodJets>=N) return true;
						      }
						   }
						}
					     }
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
      }
   return false;

}

bool OHltTree::OpenHltNPFJetPassedEleRemoval(
      int N,
      float jetPt,
      float jetEta,
      float drcut,
      float Et,
      int L1iso,
      float Tisobarrel,
      float Tisoendcap,
      float Tisoratiobarrel,
      float Tisoratioendcap,
      float HisooverETbarrel,
      float HisooverETendcap,
      float EisooverETbarrel,
      float EisooverETendcap,
      float hoverebarrel,
      float hovereendcap,
      float clusshapebarrel,
      float clusshapeendcap,
      float r9barrel,
      float r9endcap,
      float detabarrel,
      float detaendcap,
      float dphibarrel,
      float dphiendcap)
{

   float barreleta = 1.479;
   float endcapeta = 2.65;

   // Loop over all oh electrons
   for (int i=0; i<NohEle; i++)
      {
       	int isbarrel = 0;
        int isendcap = 0;
        if (TMath::Abs(ohEleEta[i]) < barreleta)
            isbarrel = 1;
        if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) < endcapeta)
               isendcap = 1;

        if (ohEleEt[i] > Et)
            {
               if (TMath::Abs(ohEleEta[i]) < endcapeta)
               {
                  if (ohEleNewSC[i]<=1)
                  {
                     if (ohElePixelSeeds[i]>0)
                     {
                      	if (ohEleL1iso[i] >= L1iso)
                        { // L1iso is 0 or 1
                           if (ohEleL1Dupl[i] == false)
                           { // remove double-counted L1 SCs
                              if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETbarrel)) ||
                                   (isendcap && ((ohEleHiso[i]/ohEleEt[i]) < HisooverETendcap)))
                              {
                                 if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETbarrel)) ||
                                    (isendcap  && ((ohEleEiso[i]/ohEleEt[i]) < EisooverETendcap)))
                                 {
                                    if ( ((isbarrel) && (ohEleHforHoverE[i]/ohEleE[i] < hoverebarrel)) ||
                                         ((isendcap) && (ohEleHforHoverE[i]/ohEleE[i] < hovereendcap)))
                                    {
                                       if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel && ohEleTiso[i] != -999.) ||
                                                           (Tisobarrel == 999.))))
                                             || (isendcap && (((ohEleTiso[i] < Tisoendcap && ohEleTiso[i] != -999.) ||
                                                               (Tisoendcap  == 999.)))))
                                       {
                                          if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i] < Tisoratiobarrel)) ||
                                              ((isendcap) && (ohEleTiso[i]/ohEleEt[i] < Tisoratioendcap)))
                                          {
                                             if ( (isbarrel && ohEleClusShap[i]  < clusshapebarrel)
                                                   || (isendcap && ohEleClusShap[i] < clusshapeendcap))
                                             {
                                                if ( (isbarrel && ohEleR9[i] < r9barrel) ||
                                                     (isendcap && ohEleR9[i] < r9endcap))
                                                {
                                                   if ( (isbarrel && TMath::Abs(ohEleDeta[i]) < detabarrel) ||
                                                        (isendcap && TMath::Abs(ohEleDeta[i]) < detaendcap))
                                                   {
                                                      if ( (isbarrel && ohEleDphi[i] < dphibarrel) ||
                                                           (isendcap && ohEleDphi[i] < dphiendcap))
                                                      {

                                                        int nGoodJets = 0;
                                                        //Loop over all oh PF jets
                                                        for (int j = 0; j < NohPFJet; j++)
                                                        {
                                                          if (pfJetPt[j] > jetPt && fabs(pfJetEta[j]) < jetEta )
                                                            { // PT, eta and JetID cuts
                                                               double deltaphi = fabs(pfJetPhi[j] -ohElePhi[i]);
                                                               if (deltaphi > 3.14159) deltaphi = (2.0 * 3.14159) - deltaphi;
                                                               double deltaRJetEle = sqrt((pfJetEta[j]-ohEleEta[i])
                                                                                   *(pfJetEta[j]-ohEleEta[i])
                                                                                   + (deltaphi *deltaphi));

                                                               if (!(deltaRJetEle < drcut)) nGoodJets++;
                                                            }
                                                        }
                                                        if (nGoodJets>=N) return true;
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
      }
   return false;

}


int OHltTree::OpenHltRPassed(
			     float Rmin,
			     float MRmin,
			     int NJmax,
			     float jetPt)
{

  //make a list of the jets
  vector<TLorentzVector*> JETS;
	
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (fabs(ohJetCorCalEta[i])>=3 || ohJetCorCalPt[i] < jetPt)  continue; // require jets with eta<3
      TLorentzVector* tmp = new TLorentzVector();
      tmp->SetPtEtaPhiE(
			ohJetCorCalPt[i],
			ohJetCorCalEta[i],
			ohJetCorCalPhi[i],
			ohJetCorCalE[i]);
		
      JETS.push_back(tmp);
    }
	
  int jetsSize = JETS.size();
	
  //Now make the hemispheres
  //for this simulation, we will used TLorentzVectors, although this is probably not
  //possible online
  if (jetsSize<2)
    return 0;
  if (NJmax!=-1 && jetsSize > NJmax)
    return 1;
  int N_comb = 1; // compute the number of combinations of jets possible
  for (int i = 0; i < jetsSize; i++)
    { // this is code is kept as close as possible
      N_comb *= 2; //to Chris' code for validation
    }
  TLorentzVector j1R,  j2R;
  double M_minR  = 9999999999.0;
  int j_count;
  for (int i=0; i<N_comb; i++)
    {
      TLorentzVector j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;
      while (j_count > 0)
	{
	  if (itemp/j_count == 1)
	    {
	      j_temp1 += *(JETS.at(count));
	    }
	  else
	    {
	      j_temp2 += *(JETS.at(count));
	    }
	  itemp -= j_count*(itemp/j_count);
	  j_count /= 2;
	  count++;
	}
      double M_temp = j_temp1.M2()+j_temp2.M2();
      if (M_temp < M_minR) 
	{
	  M_minR = M_temp;
	  j1R= j_temp1;
	  j2R= j_temp2;
	}
    }
	
  TVector3 met;
  met.SetPtEtaPhi(recoMetCal, 0, recoMetCalPhi);
	
  //CALCULATE MR
  j1R.SetPtEtaPhiM(j1R.Pt(), j1R.Eta(), j1R.Phi(), 0.0);
  j2R.SetPtEtaPhiM(j2R.Pt(), j2R.Eta(), j2R.Phi(), 0.0);
  
  if (j2R.Pt() > j1R.Pt())
    {
      TLorentzVector temp = j1R;
      j1R = j2R;
      j2R = temp;
    }
  
  //now we can calculate MTR
  double MTR = sqrt(0.5*(met.Mag()*(j1R.Pt()+j2R.Pt()) - met.Dot(j1R.Vect()+j2R.Vect())));
  
  double A = j1R.P();
  double B = j2R.P();
  double az = j1R.Pz();
  double bz = j2R.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(j1R.Px(),j1R.Py(),0.0);
  jbT.SetXYZ(j2R.Px(),j2R.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double MR = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
		   (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  
  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  MR *= mygamma; 
  
  if (MR>=MRmin && float(MTR)/float(MR)>=Rmin) return 1;
	
  return 0;
}

int OHltTree::OpenHltRMRPassed(
			     float Rmin,
			     float MRmin,
			     float RMRmin,
			     float ROffset,
			     float MROffset,
			     int NJmax,
			     float jetPt)
{

  //make a list of the jets
  vector<TLorentzVector*> JETS;
	
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (fabs(ohJetCorCalEta[i])>=3 || ohJetCorCalPt[i] < jetPt)  continue; // require jets with eta<3
      TLorentzVector* tmp = new TLorentzVector();
      tmp->SetPtEtaPhiE(
			ohJetCorCalPt[i],
			ohJetCorCalEta[i],
			ohJetCorCalPhi[i],
			ohJetCorCalE[i]);
		
      JETS.push_back(tmp);
    }
	
  int jetsSize = JETS.size();
	
  //Now make the hemispheres
  //for this simulation, we will used TLorentzVectors, although this is probably not
  //possible online
  if (jetsSize<2)
    return 0;
  if (NJmax!=-1 && jetsSize > NJmax)
    return 1;
  int N_comb = 1; // compute the number of combinations of jets possible
  for (int i = 0; i < jetsSize; i++)
    { // this is code is kept as close as possible
      N_comb *= 2; //to Chris' code for validation
    }
  TLorentzVector j1R,  j2R;
  double M_minR  = 9999999999.0;
  int j_count;
  for (int i=0; i<N_comb; i++)
    {
      TLorentzVector j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;
      while (j_count > 0)
	{
	  if (itemp/j_count == 1)
	    {
	      j_temp1 += *(JETS.at(count));
	    }
	  else
	    {
	      j_temp2 += *(JETS.at(count));
	    }
	  itemp -= j_count*(itemp/j_count);
	  j_count /= 2;
	  count++;
	}
      double M_temp = j_temp1.M2()+j_temp2.M2();
      if (M_temp < M_minR) 
	{
	  M_minR = M_temp;
	  j1R= j_temp1;
	  j2R= j_temp2;
	}
    }
	
  TVector3 met;
  met.SetPtEtaPhi(recoMetCal, 0, recoMetCalPhi);
	
  //CALCULATE MR
  j1R.SetPtEtaPhiM(j1R.Pt(), j1R.Eta(), j1R.Phi(), 0.0);
  j2R.SetPtEtaPhiM(j2R.Pt(), j2R.Eta(), j2R.Phi(), 0.0);
  
  if (j2R.Pt() > j1R.Pt())
    {
      TLorentzVector temp = j1R;
      j1R = j2R;
      j2R = temp;
    }
  
  //now we can calculate MTR
  double MTR = sqrt(0.5*(met.Mag()*(j1R.Pt()+j2R.Pt()) - met.Dot(j1R.Vect()+j2R.Vect())));
  
  double A = j1R.P();
  double B = j2R.P();
  double az = j1R.Pz();
  double bz = j2R.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(j1R.Px(),j1R.Py(),0.0);
  jbT.SetXYZ(j2R.Px(),j2R.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double MR = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
		   (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  
  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  MR *= mygamma; 
  
  float R = float(MTR)/float(MR);
  if (MR>=MRmin && R>=Rmin &&
      ((R*R-ROffset)*(MR-MROffset)>RMRmin) ) return 1;
	
  return 0;
}




int OHltTree::OpenHlt1MuonPassed(
				 double ptl1,
				 double ptl2,
				 double ptl3,
				 double dr,
				 int iso,
				 double etal2,
				 double etal3,
				 int minNHits,
				 int minNStats
				 )
{
  // This example implements the new (CMSSW_2_X) flat muon pT cuts.
  // To emulate the old behavior, the cuts should be written
  // L2:        ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]
  // L3:        ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]
	
  int rcL1 = 0;
  int rcL2 = 0;
  int rcL3 = 0;
  int rcL1L2L3 = 0;
  int NL1Mu = 8;
  int L1MinimalQuality = 1;
  int L1MaximalQuality = 7;
  int doL1L2matching = 0;
	
  for (int ic = 0; ic < 10; ic++)
    L3MuCandIDForOnia[ic] = -1;
	
  // Loop over all oh L3 muons and apply cuts
  for (int i=0; i<NohMuL3; i++)
    {
      int bestl1l2drmatchind = -1;
      double bestl1l2drmatch = 999.0;
		
      if (fabs(ohMuL3Eta[i]) < etal3)
	{ // L3 eta cut  
	  if (ohMuL3Pt[i] > ptl3)
	    { // L3 pT cut        
	      if (ohMuL3Dr[i] < dr)
		{ // L3 DR cut
		  if (ohMuL3Iso[i] >= iso)
		    { // L3 isolation
		      rcL3++;
						
		      // Begin L2 muons here. 
		      // Get best L2<->L3 match, then 
		      // begin applying cuts to L2
		      int j = ohMuL3L2idx[i]; // Get best L2<->L3 match
						
		      if ( (fabs(ohMuL2Eta[j])<etal2))
			{ // L2 eta cut
			  ////add Nhits, Nstat
			  if ((fabs(ohMuL2Eta[j])<0.9) || (fabs(ohMuL2Eta[j])>1.5 && fabs(ohMuL2Eta[j])<2.1) || ((ohMuL2Nhits[j]>= minNHits)&&(ohMuL2Nstat[j]>= minNStats)))
			    {
			    
			  if (ohMuL2Pt[j] > ptl2)
			    { // L2 pT cut
			      if (ohMuL2Iso[j] >= iso)
				{ // L2 isolation
				  rcL2++;
									
				  // Begin L1 muons here.
				  // Require there be an L1Extra muon Delta-R
				  // matched to the L2 candidate, and that it have 
				  // good quality and pass nominal L1 pT cuts 
				  for (int k = 0; k < NL1Mu; k++)
				    {
				      if ( (L1MuPt[k] < ptl1)) // L1 pT cut
					continue;
										
				      double deltaphi = fabs(ohMuL2Phi[j]-L1MuPhi[k]);
				      if (deltaphi > 3.14159)
					deltaphi = (2.0 * 3.14159) - deltaphi;
										
				      double deltarl1l2 =
					sqrt((ohMuL2Eta[j]-L1MuEta[k])
					     *(ohMuL2Eta[j]-L1MuEta[k])
					     + (deltaphi*deltaphi));
				      if (deltarl1l2 < bestl1l2drmatch)
					{
					  bestl1l2drmatchind = k;
					  bestl1l2drmatch = deltarl1l2;
					}
				    } // End loop over L1Extra muons
									
				  if (doL1L2matching == 1)
				    {
				      // Cut on L1<->L2 matching and L1 quality
				      if ((bestl1l2drmatch > 0.3)
					  || (L1MuQal[bestl1l2drmatchind]
					      < L1MinimalQuality)
					  || (L1MuQal[bestl1l2drmatchind]
					      > L1MaximalQuality))
					{
					  rcL1 = 0;
					  cout << "Failed L1-L2 match/quality" << endl;
					  cout << "L1-L2 delta-eta = "
					       << L1MuEta[bestl1l2drmatchind] << ", "
					       << ohMuL2Eta[j] << endl;
					  cout << "L1-L2 delta-pho = "
					       << L1MuPhi[bestl1l2drmatchind] << ", "
					       << ohMuL2Phi[j] << endl;
					  cout << "L1-L2 delta-R = " << bestl1l2drmatch
					       << endl;
					}
				      else
					{
					  cout << "Passed L1-L2 match/quality" << endl;
					  L3MuCandIDForOnia[rcL1L2L3] = i;
					  rcL1++;
					  rcL1L2L3++;
					} // End L1 matching and quality cuts	      
				    }
				  else
				    {
				      L3MuCandIDForOnia[rcL1L2L3] = i;
				      rcL1L2L3++;
				    }
				} // End L2 isolation cut 
			    }//end Nhits, NStat cuts
			    } // End L2 eta cut
			} // End L2 pT cut
		    } // End L3 isolation cut
		} // End L3 DR cut
	    } // End L3 pT cut
	} // End L3 eta cut
    } // End loop over L3 muons		      
	
  return rcL1L2L3;
}

int OHltTree::OpenHlt1MuonPassed(
				 vector<double> muThresholds,
				 double dr,
				 int iso,
				 double etal2,
				 double etal3,
				 int minNHits,
				 int minNStats 
				 )
{
  // This example implements the new (CMSSW_2_X) flat muon pT cuts.
  // To emulate the old behavior, the cuts should be written
  // L2:        ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]
  // L3:        ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]
	
  int rcL1 = 0;
  int rcL2 = 0;
  int rcL3 = 0;
  int rcL1L2L3 = 0;
  int NL1Mu = 8;
  int L1MinimalQuality = 1;
  int L1MaximalQuality = 7;
  int doL1L2matching = 0;
	
  for (int ic = 0; ic < 10; ic++)
    L3MuCandIDForOnia[ic] = -1;
	
  // Loop over all oh L3 muons and apply cuts
  for (int i=0; i<NohMuL3; i++)
    {
      int bestl1l2drmatchind = -1;
      double bestl1l2drmatch = 999.0;
		
      if (fabs(ohMuL3Eta[i]) < etal3)
	{ // L3 eta cut  
	  if (ohMuL3Pt[i] > muThresholds[2])
	    { // L3 pT cut        
	      if (ohMuL3Dr[i] < dr)
		{ // L3 DR cut
		  if (ohMuL3Iso[i] >= iso)
		    { // L3 isolation
		      rcL3++;
						
		      // Begin L2 muons here. 
		      // Get best L2<->L3 match, then 
		      // begin applying cuts to L2
		      int j = ohMuL3L2idx[i]; // Get best L2<->L3 match
						
		      if ( (fabs(ohMuL2Eta[j])<etal2))
			{ // L2 eta cut
			  if ((fabs(ohMuL2Eta[j])<0.9) || (fabs(ohMuL2Eta[j])>1.5 && fabs(ohMuL2Eta[j])<2.1) || ((ohMuL2Nhits[j]>= minNHits)&&(ohMuL2Nstat[j]>minNStats)))
			    {
			  if (ohMuL2Pt[j] >  muThresholds[1] )
			    { // L2 pT cut
			      if (ohMuL2Iso[j] >= iso)
				{ // L2 isolation
				  rcL2++;
									
				  // Begin L1 muons here.
				  // Require there be an L1Extra muon Delta-R
				  // matched to the L2 candidate, and that it have 
				  // good quality and pass nominal L1 pT cuts 
				  for (int k = 0; k < NL1Mu; k++)
				    {
				      if ( (L1MuPt[k] < muThresholds[0])) // L1 pT cut
					continue;
										
				      double deltaphi = fabs(ohMuL2Phi[j]-L1MuPhi[k]);
				      if (deltaphi > 3.14159)
					deltaphi = (2.0 * 3.14159) - deltaphi;
										
				      double deltarl1l2 =
					sqrt((ohMuL2Eta[j]-L1MuEta[k])
					     *(ohMuL2Eta[j]-L1MuEta[k])
					     + (deltaphi*deltaphi));
				      if (deltarl1l2 < bestl1l2drmatch)
					{
					  bestl1l2drmatchind = k;
					  bestl1l2drmatch = deltarl1l2;
					}
				    } // End loop over L1Extra muons
									
				  if (doL1L2matching == 1)
				    {
				      // Cut on L1<->L2 matching and L1 quality
				      if ((bestl1l2drmatch > 0.3)
					  || (L1MuQal[bestl1l2drmatchind]
					      < L1MinimalQuality)
					  || (L1MuQal[bestl1l2drmatchind]
					      > L1MaximalQuality))
					{
					  rcL1 = 0;
					  cout << "Failed L1-L2 match/quality" << endl;
					  cout << "L1-L2 delta-eta = "
					       << L1MuEta[bestl1l2drmatchind] << ", "
					       << ohMuL2Eta[j] << endl;
					  cout << "L1-L2 delta-pho = "
					       << L1MuPhi[bestl1l2drmatchind] << ", "
					       << ohMuL2Phi[j] << endl;
					  cout << "L1-L2 delta-R = " << bestl1l2drmatch
					       << endl;
					}
				      else
					{
					  cout << "Passed L1-L2 match/quality" << endl;
					  L3MuCandIDForOnia[rcL1L2L3] = i;
					  rcL1++;
					  rcL1L2L3++;
					} // End L1 matching and quality cuts	      
				    }
				  else
				    {
				      L3MuCandIDForOnia[rcL1L2L3] = i;
				      rcL1L2L3++;
				    }
				} // End L2 isolation cut 
			    }//end NHits, NStat
			    } // End L2 eta cut
			} // End L2 pT cut
		    } // End L3 isolation cut
		} // End L3 DR cut
	    } // End L3 pT cut
	} // End L3 eta cut
    } // End loop over L3 muons		      
	
  return rcL1L2L3;
}


//// Separating between Jets and Muon.....
int OHltTree::OpenHlt1MuonIsoJetPassed(
				       double ptl1,
				       double ptl2,
				       double ptl3,
				       double dr,
				       int iso,
				       double JetPt,
				       double JetEta)
{
  // This example implements the new (CMSSW_2_X) flat muon pT cuts.
  // To emulate the old behavior, the cuts should be written
  // L2:        ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]
  // L3:        ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]
	
  int rcL1 = 0;
  int rcL2 = 0;
  int rcL3 = 0;
  int rcL1L2L3 = 0;
  int NL1Mu = 8;
  int L1MinimalQuality = 4;
  int L1MaximalQuality = 7;
  int doL1L2matching = 0;
	
  for (int ic = 0; ic < 10; ic++)
    L3MuCandIDForOnia[ic] = -1;
	
  // Loop over all oh L3 muons and apply cuts
  for (int i=0; i<NohMuL3; i++)
    {
      int bestl1l2drmatchind = -1;
      double bestl1l2drmatch = 999.0;
		
      if (fabs(ohMuL3Eta[i]) < 2.5)
	{ // L3 eta cut  
	  if (ohMuL3Pt[i] > ptl3)
	    { // L3 pT cut        
	      if (ohMuL3Dr[i] < dr)
		{ // L3 DR cut
		  if (ohMuL3Iso[i] >= iso)
		    { // L3 isolation
						
		      // Loop over all oh corrected jets    
		      float minDR = 100.;
		      for (int j=0; j <NohJetCorCal; j++)
			{
			  if (ohJetCorCalPt[j]>JetPt && fabs(ohJetCorCalEta[j])
			      <JetEta)
			    { // Jet pT cut
			      double deltaphi =
				fabs(ohJetCorCalPhi[j]-ohMuL3Phi[i]);
			      if (deltaphi > 3.14159)
				deltaphi = (2.0 * 3.14159) - deltaphi;
			      float deltaRMuJet = sqrt((recoJetCorCalEta[j]
							-ohMuL3Eta[i])*(recoJetCorCalEta[j]-ohMuL3Eta[i])
						       + (deltaphi*deltaphi));
			      if (deltaRMuJet < minDR)
				{
				  minDR = deltaRMuJet;
				}
			    }
			}
		      if (minDR < 0.3)
			break;
						
		      rcL3++;
						
		      // Begin L2 muons here. 
		      // Get best L2<->L3 match, then 
		      // begin applying cuts to L2
		      int j = ohMuL3L2idx[i]; // Get best L2<->L3 match
						
		      if ( (fabs(ohMuL2Eta[j])<2.5))
			{ // L2 eta cut
			  if (ohMuL2Pt[j] > ptl2)
			    { // L2 pT cut
			      if (ohMuL2Iso[j] >= iso)
				{ // L2 isolation
				  rcL2++;
									
				  // Begin L1 muons here.
				  // Require there be an L1Extra muon Delta-R
				  // matched to the L2 candidate, and that it have 
				  // good quality and pass nominal L1 pT cuts 
				  for (int k = 0; k < NL1Mu; k++)
				    {
				      if ( (L1MuPt[k] < ptl1)) // L1 pT cut
					continue;
										
				      double deltaphi = fabs(ohMuL2Phi[j]-L1MuPhi[k]);
				      if (deltaphi > 3.14159)
					deltaphi = (2.0 * 3.14159) - deltaphi;
										
				      double deltarl1l2 =
					sqrt((ohMuL2Eta[j]-L1MuEta[k])
					     *(ohMuL2Eta[j]-L1MuEta[k])
					     + (deltaphi*deltaphi));
				      if (deltarl1l2 < bestl1l2drmatch)
					{
					  bestl1l2drmatchind = k;
					  bestl1l2drmatch = deltarl1l2;
					}
				    } // End loop over L1Extra muons
									
				  if (doL1L2matching == 1)
				    {
				      // Cut on L1<->L2 matching and L1 quality
				      if ((bestl1l2drmatch > 0.3)
					  || (L1MuQal[bestl1l2drmatchind]
					      < L1MinimalQuality)
					  || (L1MuQal[bestl1l2drmatchind]
					      > L1MaximalQuality))
					{
					  rcL1 = 0;
					  cout << "Failed L1-L2 match/quality" << endl;
					  cout << "L1-L2 delta-eta = "
					       << L1MuEta[bestl1l2drmatchind] << ", "
					       << ohMuL2Eta[j] << endl;
					  cout << "L1-L2 delta-pho = "
					       << L1MuPhi[bestl1l2drmatchind] << ", "
					       << ohMuL2Phi[j] << endl;
					  cout << "L1-L2 delta-R = " << bestl1l2drmatch
					       << endl;
					}
				      else
					{
					  cout << "Passed L1-L2 match/quality" << endl;
					  L3MuCandIDForOnia[rcL1L2L3] = i;
					  rcL1++;
					  rcL1L2L3++;
					} // End L1 matching and quality cuts      
				    }
				  else
				    {
				      L3MuCandIDForOnia[rcL1L2L3] = i;
				      rcL1L2L3++;
				    }
				} // End L2 isolation cut 
			    } // End L2 eta cut
			} // End L2 pT cut
		    } // End L3 isolation cut
		} // End L3 DR cut
	    } // End L3 pT cut
	} // End L3 eta cut
    } // End loop over L3 muons      
  return rcL1L2L3;
	
}

int OHltTree::OpenHlt1L2MuonPassed(double ptl1, double ptl2, double dr)
{
  // This is a modification of the standard Hlt1Muon code, which does not consider L3 information 
	
  int rcL1 = 0;
  int rcL2 = 0;
  int rcL1L2L3 = 0;
  int NL1Mu = 8;
  int L1MinimalQuality = 3;
  int L1MaximalQuality = 7;
  int doL1L2matching = 0;
	
  // Loop over all oh L2 muons and apply cuts 
  for (int j=0; j<NohMuL2; j++)
    {
      int bestl1l2drmatchind = -1;
      double bestl1l2drmatch = 999.0;
		
      if (fabs(ohMuL2Eta[j])>=2.5)
	continue; // L2 eta cut 
      if (ohMuL2Pt[j] <= ptl2)
	continue; // L2 pT cut 
      rcL2++;
		
      // Begin L1 muons here. 
      // Require there be an L1Extra muon Delta-R 
      // matched to the L2 candidate, and that it have  
      // good quality and pass nominal L1 pT cuts  
      for (int k = 0; k < NL1Mu; k++)
	{
	  if ( (L1MuPt[k] < ptl1))
	    continue; // L1 pT cut 
			
	  double deltaphi = fabs(ohMuL2Phi[j]-L1MuPhi[k]);
	  if (deltaphi > 3.14159)
	    deltaphi = (2.0 * 3.14159) - deltaphi;
	  double deltarl1l2 = sqrt((ohMuL2Eta[j]-L1MuEta[k])*(ohMuL2Eta[j]
							      -L1MuEta[k]) + (deltaphi*deltaphi));
	  if (deltarl1l2 < bestl1l2drmatch)
	    {
	      bestl1l2drmatchind = k;
	      bestl1l2drmatch = deltarl1l2;
	    }
	} // End loop over L1Extra muons 
      if (doL1L2matching == 1)
	{
	  // Cut on L1<->L2 matching and L1 quality 
	  if ((bestl1l2drmatch > 0.3) || (L1MuQal[bestl1l2drmatchind]
					  < L1MinimalQuality) || (L1MuQal[bestl1l2drmatchind]
								  > L1MaximalQuality))
	    {
	      rcL1 = 0;
	      cout << "Failed L1-L2 match/quality" << endl;
	      cout << "L1-L2 delta-eta = " << L1MuEta[bestl1l2drmatchind] << ", "
		   << ohMuL2Eta[j] << endl;
	      cout << "L1-L2 delta-pho = " << L1MuPhi[bestl1l2drmatchind] << ", "
		   << ohMuL2Phi[j] << endl;
	      cout << "L1-L2 delta-R = " << bestl1l2drmatch << endl;
	    }
	  else
	    {
	      cout << "Passed L1-L2 match/quality" << endl;
	      rcL1++;
	      rcL1L2L3++;
	    } // End L1 matching and quality cuts            
	}
      else
	{
	  rcL1L2L3++;
	}
    } // End L2 loop over muons 
  return rcL1L2L3;
}

int OHltTree::OpenHlt1L2MuonNoVertexPassed(double ptl1, double ptl2, double dr) 
{ 
  // This is a modification of the standard Hlt1Muon code, which does not consider L3 information  
	
  int rcL1 = 0; 
  int rcL2 = 0; 
  int rcL1L2L3 = 0; 
  int NL1Mu = 8; 
  int L1MinimalQuality = 3; 
  int L1MaximalQuality = 7; 
  int doL1L2NoVtxmatching = 0; 
	
  // Loop over all oh L2NoVtx muons and apply cuts  
  for (int j=0; j<NohMuL2NoVtx; j++) 
    { 
      int bestl1l2drmatchind = -1; 
      double bestl1l2drmatch = 999.0; 
		
      if (fabs(ohMuL2NoVtxEta[j])>=2.5) 
	continue; // L2NoVtx eta cut  
      if (ohMuL2NoVtxPt[j] <= ptl2) 
	continue; // L2NoVtx pT cut  
      rcL2++; 
		
      // Begin L1 muons here.  
      // Require there be an L1Extra muon Delta-R  
      // matched to the L2NoVtx candidate, and that it have   
      // good quality and pass nominal L1 pT cuts   
      for (int k = 0; k < NL1Mu; k++) 
	{ 
	  if ( (L1MuPt[k] < ptl1)) 
	    continue; // L1 pT cut  
			
	  double deltaphi = fabs(ohMuL2NoVtxPhi[j]-L1MuPhi[k]); 
	  if (deltaphi > 3.14159) 
	    deltaphi = (2.0 * 3.14159) - deltaphi; 
	  double deltarl1l2 = sqrt((ohMuL2NoVtxEta[j]-L1MuEta[k])*(ohMuL2NoVtxEta[j] 
								   -L1MuEta[k]) + (deltaphi*deltaphi)); 
	  if (deltarl1l2 < bestl1l2drmatch) 
	    { 
	      bestl1l2drmatchind = k; 
	      bestl1l2drmatch = deltarl1l2; 
	    } 
	} // End loop over L1Extra muons  
      if (doL1L2NoVtxmatching == 1) 
	{ 
	  // Cut on L1<->L2NoVtx matching and L1 quality  
	  if ((bestl1l2drmatch > 0.3) || (L1MuQal[bestl1l2drmatchind] 
					  < L1MinimalQuality) || (L1MuQal[bestl1l2drmatchind] 
								  > L1MaximalQuality)) 
	    { 
	      rcL1 = 0; 
	    } 
	  else 
	    { 
	      cout << "Passed L1-L2NoVtx match/quality" << endl; 
	      rcL1++; 
	      rcL1L2L3++; 
	    } // End L1 matching and quality cuts             
	} 
      else 
	{ 
	  rcL1L2L3++; 
	} 
    } // End L2 loop over muons  
  return rcL1L2L3; 
} 


int OHltTree::OpenHltMuPixelPassed(
				   double ptPix,
				   double pPix,
				   double etaPix,
				   double DxyPix,
				   double DzPix,
				   int NHitsPix,
				   double normChi2Pix,
				   double *massMinPix,
				   double *massMaxPix,
				   double DzMuonPix,
				   bool checkChargePix)
{
	
  //   printf("\n\n");
  const double muMass = 0.105658367;
  TLorentzVector pix4Mom, mu4Mom, onia4Mom;
  int iNPix = 0;
  //reset counter variables:
  for (int iMu = 0; iMu < 10; iMu++)
    {
      L3PixelCandIDForOnia[iMu] = -1;
      L3MuPixCandIDForOnia[iMu] = -1;
    }
	
  //0.) check how many L3 muons there are:
  int nMuons = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuCandIDForOnia[iMu] > -1)
      nMuons++;
	
  //1.) loop over the Pixel tracks
  for (int iP = 0; iP < NohOniaPixel; iP++)
    {
		
      //select those that survive the kinematical and
      //topological selection cuts
      if (fabs(ohOniaPixelEta[iP]) > etaPix)
	continue; //eta cut
      if (ohOniaPixelPt[iP] < ptPix)
	continue; //pT cut
      double momThisPix = ohOniaPixelPt[iP] * cosh(ohOniaPixelEta[iP]);
      if (momThisPix < pPix)
	continue; //momentum cut
      if (ohOniaPixelHits[iP] < NHitsPix)
	continue; //min. nb. of hits
      if (ohOniaPixelNormChi2[iP] > normChi2Pix)
	continue; //chi2 cut
      if (fabs(ohOniaPixelDr[iP]) > DxyPix)
	continue; //Dr cut
      if (fabs(ohOniaPixelDz[iP]) > DzPix)
	continue;
		
      pix4Mom.SetPtEtaPhiM(
			   ohOniaPixelPt[iP],
			   ohOniaPixelEta[iP],
			   ohOniaPixelPhi[iP],
			   muMass);
      //2.) loop now over all L3 muons and check if they would give a
      //Onia (J/psi or upsilon) pair:
      for (int iMu = 0; iMu < nMuons; iMu++)
	{
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuCandIDForOnia[iMu]],
			      ohMuL3Eta[L3MuCandIDForOnia[iMu]],
			      ohMuL3Phi[L3MuCandIDForOnia[iMu]],
			      muMass);
	  onia4Mom = pix4Mom + mu4Mom;
			
	  double oniaMass = onia4Mom.M();
	  if (oniaMass < massMinPix[0] || oniaMass> massMaxPix[1]) continue; //mass cut
	  if(oniaMass> massMaxPix[0] && oniaMass < massMinPix[1]) continue; //mass cut
	  if(checkChargePix)
	    if(ohMuL3Chg[iMu] == ohOniaPixelChg[iP]) continue; //charge cut
	  if(fabs(ohMuL3Dz[iMu] - ohOniaPixelDz[iP])> DzMuonPix) continue;
			
	  //store the surviving pixel-muon combinations:
	  if(iNPix < 10)
	    {
	      L3PixelCandIDForOnia[iNPix] = iP;
	      L3MuPixCandIDForOnia[iNPix] = iMu;
	      iNPix++;
	    }
	  //       printf("mu[%d]-pixel[%d] inv. mass %f\n",
	  //           L3MuCandIDForOnia[iMu], iP, oniaMass);
	}
    }
	
  //   hNPixelCand->Fill(iNPix);
  return iNPix;
}

int OHltTree::OpenHltMuTrackPassed(
				   double ptTrack,
				   double pTrack,
				   double etaTrack,
				   double DxyTrack,
				   double DzTrack,
				   int NHitsTrack,
				   double normChi2Track,
				   double *massMinTrack,
				   double *massMaxTrack,
				   double DzMuonTrack,
				   bool checkChargeTrack)
{
	
  double pixMatchingDeltaR = 0.03;
  const double muMass = 0.105658367;
  TLorentzVector track4Mom, mu4Mom, onia4Mom;
  int iNTrack = 0;
	
  //0.) check how many pixel-muon combinations there are:
  int nComb = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuPixCandIDForOnia[iMu] > -1)
      nComb++;
	
  //   printf("OpenHltMuTrackPassed: %d incoming pixels and %d tracks\n", nComb, NohOniaTrack);
	
  //1.) loop over the Tracker tracks
  for (int iT = 0; iT < NohOniaTrack; iT++)
    {
		
      //select those that survive the kinematical and
      //topological selection cuts
      if (fabs(ohOniaTrackEta[iT]) > etaTrack)
	continue; //eta cut
      if (ohOniaTrackPt[iT] < ptTrack)
	continue; //pT cut
      double momThisTrack = ohOniaTrackPt[iT] * cosh(ohOniaTrackEta[iT]);
      //     printf("track[%d] has eta %f, pT %f and mom %f\n",
      //         iT, ohOniaTrackEta[iT], ohOniaTrackPt[iT], momThisTrack);
      if (momThisTrack < pTrack)
	continue; //momentum cut
      if (ohOniaTrackHits[iT] < NHitsTrack)
	continue; //min. nb. of hits
      if (ohOniaTrackNormChi2[iT] > normChi2Track)
	continue; //chi2 cut
      if (fabs(ohOniaTrackDr[iT]) > DxyTrack)
	continue; //Dr cut
      if (fabs(ohOniaTrackDz[iT]) > DzTrack)
	continue;
		
      //2.) loop over the pixels candidates to see whether the track
      //under investigation has a match to the pixel track
      bool trackMatched = false;
      for (int iPix = 0; iPix < nComb; iPix++)
	{
			
	  if (trackMatched)
	    break; //in case the track was already matched
	  if (L3PixelCandIDForOnia[iPix] < 0)
	    continue; //in case the pixel has been matched to a previous track
			
	  double deltaEta = ohOniaPixelEta[L3PixelCandIDForOnia[iPix]]
	    - ohOniaTrackEta[iT];
	  double deltaPhi = ohOniaPixelPhi[L3PixelCandIDForOnia[iPix]]
	    - ohOniaTrackPhi[iT];
	  double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
			
	  if (deltaR > pixMatchingDeltaR)
	    continue;
	  //       printf("track[%d], pixel[%d], delta R %f\n", iT, L3PixelCandIDForOnia[iPix], deltaR);
			
	  trackMatched = true;
	  L3PixelCandIDForOnia[iPix] = -1; //deactivate this candidate to not match it to any further track
			
	  track4Mom.SetPtEtaPhiM(
				 ohOniaTrackPt[iT],
				 ohOniaTrackEta[iT],
				 ohOniaTrackPhi[iT],
				 muMass);
	  //check if the matched tracker track combined with the
	  //muon gives again an opposite sign onia:
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuPixCandIDForOnia[iPix]],
			      ohMuL3Eta[L3MuPixCandIDForOnia[iPix]],
			      ohMuL3Phi[L3MuPixCandIDForOnia[iPix]],
			      muMass);
	  onia4Mom = track4Mom + mu4Mom;
			
	  double oniaMass = onia4Mom.M();
	  //       printf("mu[%d]-track[%d] inv. mass %f\n",
	  //           L3MuPixCandIDForOnia[iPix], iT, oniaMass);
			
	  if (oniaMass < massMinTrack[0] || oniaMass> massMaxTrack[1]) continue; //mass cut
	  if(oniaMass> massMaxTrack[0] && oniaMass < massMinTrack[1]) continue; //mass cut
			
	  //       printf("surviving: mu[%d]-track[%d] inv. mass %f\n",
	  //           L3MuPixCandIDForOnia[iPix], iT, oniaMass);
			
	  if(checkChargeTrack)
	    if(ohMuL3Chg[L3MuPixCandIDForOnia[iPix]] == ohOniaTrackChg[iT]) continue; //charge cut
	  if(fabs(ohMuL3Dz[L3MuPixCandIDForOnia[iPix]] - ohOniaTrackDz[iT])> DzMuonTrack) continue; //deltaZ cut
			
	  //store the surviving track-muon combinations:
	  if(iNTrack < 10)
	    iNTrack++;
			
	  break; //don't check further pixels... go to next track
	}
    }
	
  //   if(iNTrack > 0)
  //     printf("found %d final candidates!!!\n", iNTrack);
  return iNTrack;
}

int OHltTree::OpenHltMuPixelPassed_JPsi(
					double ptPix,
					double pPix,
					double etaPix,
					double DxyPix,
					double DzPix,
					int NHitsPix,
					double normChi2Pix,
					double *massMinPix,
					double *massMaxPix,
					double DzMuonPix,
					bool checkChargePix,
					int histIndex)
{
  //   printf("in OpenHltMuPixelPassed_JPsi \n\n");
  const double muMass = 0.105658367;
  TLorentzVector pix4Mom, mu4Mom, onia4Mom;
  int iNPix = 0;
  //reset counter variables:
  for (int iMu = 0; iMu < 10; iMu++)
    {
      L3PixelCandIDForJPsi[iMu] = -1;
      L3MuPixCandIDForJPsi[iMu] = -1;
    }
	
  //0.) check how many L3 muons there are:
  int nMuons = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuCandIDForOnia[iMu] > -1)
      nMuons++;
	
  //   Int_t countCut = 0, countOniaCut = 0;
  //1.) loop over the Pixel tracks
  for (int iP = 0; iP < NohOniaPixel; iP++)
    {
		
      //     countCut = 0;
      //     hEta[histIndex][0][countCut]->Fill(ohOniaPixelEta[iP]);
      //     hPt[histIndex][0][countCut]->Fill(ohOniaPixelPt[iP]);
      //     hHits[histIndex][0][countCut]->Fill(ohOniaPixelHits[iP]);
      //     hNormChi2[histIndex][0][countCut]->Fill(ohOniaPixelNormChi2[iP]);
      //     hDxy[histIndex][0][countCut]->Fill(ohOniaPixelDr[iP]);
      //     hDz[histIndex][0][countCut]->Fill(ohOniaPixelDz[iP]);
		
      //select those that survive the kinematical and
      //topological selection cuts
      if (fabs(ohOniaPixelEta[iP]) > etaPix)
	continue; //eta cut
      if (ohOniaPixelPt[iP] < ptPix)
	continue; //pT cut
		
      double momThisPix = ohOniaPixelPt[iP] * cosh(ohOniaPixelEta[iP]);
      //     hP[histIndex][0][countCut]->Fill(momThisPix);
      //     countCut++;
		
      if (momThisPix < pPix)
	continue; //momentum cut
      if (ohOniaPixelHits[iP] < NHitsPix)
	continue; //min. nb. of hits
      if (ohOniaPixelNormChi2[iP] > normChi2Pix)
	continue; //chi2 cut
      if (fabs(ohOniaPixelDr[iP]) > DxyPix)
	continue; //Dr cut
      if (fabs(ohOniaPixelDz[iP]) > DzPix)
	continue;
		
      //     hEta[histIndex][0][countCut]->Fill(ohOniaPixelEta[iP]);
      //     hPt[histIndex][0][countCut]->Fill(ohOniaPixelPt[iP]);
      //     hHits[histIndex][0][countCut]->Fill(ohOniaPixelHits[iP]);
      //     hNormChi2[histIndex][0][countCut]->Fill(ohOniaPixelNormChi2[iP]);
      //     hDxy[histIndex][0][countCut]->Fill(ohOniaPixelDr[iP]);
      //     hDz[histIndex][0][countCut]->Fill(ohOniaPixelDz[iP]);
      //     hP[histIndex][0][countCut]->Fill(momThisPix);
      //     countCut++;
		
      pix4Mom.SetPtEtaPhiM(
			   ohOniaPixelPt[iP],
			   ohOniaPixelEta[iP],
			   ohOniaPixelPhi[iP],
			   muMass);
      //2.) loop now over all L3 muons and check if they would give a
      //JPsi pair:
      for (int iMu = 0; iMu < nMuons; iMu++)
	{
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuCandIDForOnia[iMu]],
			      ohMuL3Eta[L3MuCandIDForOnia[iMu]],
			      ohMuL3Phi[L3MuCandIDForOnia[iMu]],
			      muMass);
	  onia4Mom = pix4Mom + mu4Mom;
			
	  double oniaMass = onia4Mom.M();
	  //       printf("mu[%d]-pixel[%d] inv. mass %f\n",
	  //           L3MuCandIDForOnia[iMu], iP, oniaMass);
	  //       countOniaCut = 0;
	  //       if(oniaMass > 5.0) continue; //Only JPsi 
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if (oniaMass < massMinPix[0] || oniaMass> massMaxPix[0]) continue; //mass cut
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(checkChargePix)
	    if(ohMuL3Chg[iMu] == ohOniaPixelChg[iP]) continue; //charge cut
			
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(fabs(ohMuL3Dz[iMu] - ohOniaPixelDz[iP])> DzMuonPix) continue;
			
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
			
	  //store the surviving pixel-muon combinations:
	  if(iNPix < 10)
	    {
	      L3PixelCandIDForJPsi[iNPix] = iP;
	      L3MuPixCandIDForJPsi[iNPix] = iMu;
	      iNPix++;
	    }
	  //       printf("surviving: mu[%d]-pixel[%d] inv. mass %f\n",
	  //           L3MuCandIDForOnia[iMu], iP, oniaMass);
	}
    }
	
  //   hNCand[histIndex][0]->Fill(iNPix);
	
  //Pixel Eta, Pt, P, DR
  if(iNPix!=0)
    {
      for(int inP=0;inP<iNPix;inP++)
	{
			
	  //       hPixCandEta[histIndex]->Fill(ohOniaPixelEta[L3PixelCandIDForJPsi[inP]]);
	  //       hPixCandPt[histIndex]->Fill(ohOniaPixelPt[L3PixelCandIDForJPsi[inP]]);
	  //       hPixCandP[histIndex]->Fill(ohOniaPixelPt[L3PixelCandIDForJPsi[inP]] * cosh(ohOniaPixelEta[L3PixelCandIDForJPsi[inP]]));
			
	  //       if(iNPix>=2){
	  //         for(int jnP=inP+1;jnP<iNPix;jnP++){
	  //           if(inP!=jnP){
	  //              double dEta = fabs(ohOniaPixelEta[L3PixelCandIDForJPsi[inP]]-ohOniaPixelEta[L3PixelCandIDForJPsi[jnP]]);
	  //              double dPhi = fabs(ohOniaPixelPhi[L3PixelCandIDForJPsi[inP]]-ohOniaPixelPhi[L3PixelCandIDForJPsi[jnP]]);
	  //              if(dPhi>TMath::Pi()) dPhi = 2.0*TMath::Pi()-dPhi;
	  //              hPixCanddr[histIndex]->Fill(sqrt(pow(dEta,2)+pow(dPhi,2)));
	  //           }
	  //         }
	  //       }
	}
    }
	
  return iNPix;
}

int OHltTree::OpenHltMuTrackPassed_JPsi(
					double ptTrack,
					double pTrack,
					double etaTrack,
					double DxyTrack,
					double DzTrack,
					int NHitsTrack,
					double normChi2Track,
					double *massMinTrack,
					double *massMaxTrack,
					double DzMuonTrack,
					bool checkChargeTrack,
					int histIndex)
{
	
  double pixMatchingDeltaR = 0.01;
  const double muMass = 0.105658367;
  TLorentzVector track4Mom, mu4Mom, onia4Mom;
  int iNTrack = 0;
	
  //0.) check how many pixel-muon combinations there are:
  int nComb = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuPixCandIDForJPsi[iMu] > -1)
      nComb++;
	
  //   printf("OpenHltMuTrackPassed_JPsi: %d incoming pixels and %d tracks\n", nComb, NohOniaTrack);
  //   Int_t countCut = 0, countOniaCut = 0;
  //1.) loop over the Tracker tracks
  for (int iT = 0; iT < NohOniaTrack; iT++)
    {
		
      //select those that survive the kinematical and
      //topological selection cuts
      //     countCut = 0;
      //     hEta[histIndex][1][countCut]->Fill(ohOniaTrackEta[iT]);
      //     hPt[histIndex][1][countCut]->Fill(ohOniaTrackPt[iT]);
      //     hHits[histIndex][1][countCut]->Fill(ohOniaTrackHits[iT]);
      //     hNormChi2[histIndex][1][countCut]->Fill(ohOniaTrackNormChi2[iT]);
      //     hDxy[histIndex][1][countCut]->Fill(ohOniaTrackDr[iT]);
      //     hDz[histIndex][1][countCut]->Fill(ohOniaTrackDz[iT]);
		
      if (fabs(ohOniaTrackEta[iT]) > etaTrack)
	continue; //eta cut
      if (ohOniaTrackPt[iT] < ptTrack)
	continue; //pT cut
      double momThisTrack = ohOniaTrackPt[iT] * cosh(ohOniaTrackEta[iT]);
      //     printf("track[%d] has eta %f, pT %f and mom %f\n",
      //         iT, ohOniaTrackEta[iT], ohOniaTrackPt[iT], momThisTrack);
      //     hP[histIndex][1][countCut]->Fill(momThisTrack);
      //     countCut++;
		
      if (momThisTrack < pTrack)
	continue; //momentum cut
      if (ohOniaTrackHits[iT] < NHitsTrack)
	continue; //min. nb. of hits
      if (ohOniaTrackNormChi2[iT] > normChi2Track)
	continue; //chi2 cut
      if (fabs(ohOniaTrackDr[iT]) > DxyTrack)
	continue; //Dr cut
      if (fabs(ohOniaTrackDz[iT]) > DzTrack)
	continue;
		
      //     hEta[histIndex][1][countCut]->Fill(ohOniaTrackEta[iT]);
      //     hPt[histIndex][1][countCut]->Fill(ohOniaTrackPt[iT]);
      //     hHits[histIndex][1][countCut]->Fill(ohOniaTrackHits[iT]);
      //     hNormChi2[histIndex][1][countCut]->Fill(ohOniaTrackNormChi2[iT]);
      //     hDxy[histIndex][1][countCut]->Fill(ohOniaTrackDr[iT]);
      //     hDz[histIndex][1][countCut]->Fill(ohOniaTrackDz[iT]);
      //     hP[histIndex][1][countCut]->Fill(momThisTrack);
      //     countCut++;
		
      //     printf("track %d surviving kinematical pre-selection\n", iT);
      //2.) loop over the pixels candidates to see whether the track
      //under investigation has a match to the pixel track
      bool trackMatched = false;
      for (int iPix = 0; iPix < nComb; iPix++)
	{
	  if (trackMatched)
	    break; //in case the track was already matched
	  if (L3PixelCandIDForJPsi[iPix] < 0)
	    continue; //in case the pixel has been matched to a previous track
			
	  double deltaEta = ohOniaPixelEta[L3PixelCandIDForJPsi[iPix]]
	    - ohOniaTrackEta[iT];
	  double deltaPhi = ohOniaPixelPhi[L3PixelCandIDForJPsi[iPix]]
	    - ohOniaTrackPhi[iT];
	  if (deltaPhi>TMath::Pi())
	    deltaPhi = 2.0*TMath::Pi()-deltaPhi;
	  double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
			
	  //       printf("delta R = %f\n", deltaR);
	  if (deltaR > pixMatchingDeltaR)
	    continue;
	  //       printf("track[%d] and pixel[%d] are compatible (deltaR %f)\n", iT, L3PixelCandIDForJPsi[iPix], deltaR);
			
	  trackMatched = true;
	  L3PixelCandIDForJPsi[iPix] = -1; //deactivate this candidate to not match it to any further track
			
	  track4Mom.SetPtEtaPhiM(
				 ohOniaTrackPt[iT],
				 ohOniaTrackEta[iT],
				 ohOniaTrackPhi[iT],
				 muMass);
	  //check if the matched tracker track combined with the
	  //muon gives again an opposite sign onia:
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuPixCandIDForJPsi[iPix]],
			      ohMuL3Eta[L3MuPixCandIDForJPsi[iPix]],
			      ohMuL3Phi[L3MuPixCandIDForJPsi[iPix]],
			      muMass);
	  onia4Mom = track4Mom + mu4Mom;
			
	  double oniaMass = onia4Mom.M();
	  //       printf("mu[%d]-track[%d] inv. mass %f\n",
	  //           L3MuPixCandIDForJPsi[iPix], iT, oniaMass);
			
	  //       countOniaCut = 0;
	  //       if(oniaMass>5.0) continue; //Only JPsi
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if (oniaMass < massMinTrack[0] || oniaMass> massMaxTrack[0]) continue; //mass cut
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(checkChargeTrack)
	    if(ohMuL3Chg[L3MuPixCandIDForJPsi[iPix]] == ohOniaTrackChg[iT]) continue; //charge cut
			
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(fabs(ohMuL3Dz[L3MuPixCandIDForJPsi[iPix]] - ohOniaTrackDz[iT])> DzMuonTrack) continue; //deltaZ cut
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
			
	  //store the surviving track-muon combinations:
	  if(iNTrack < 10)
	    iNTrack++;
	  break; //don't check further pixels... go to next track
	}
    }
	
  //   hNCand[histIndex][1]->Fill(iNTrack);
  return iNTrack;
}

int OHltTree::OpenHltMuPixelPassed_Ups(
				       double ptPix,
				       double pPix,
				       double etaPix,
				       double DxyPix,
				       double DzPix,
				       int NHitsPix,
				       double normChi2Pix,
				       double *massMinPix,
				       double *massMaxPix,
				       double DzMuonPix,
				       bool checkChargePix,
				       int histIndex)
{
	
  const double muMass = 0.105658367;
  TLorentzVector pix4Mom, mu4Mom, onia4Mom;
  int iNPix = 0;
  //reset counter variables:
  for (int iMu = 0; iMu < 10; iMu++)
    {
      L3PixelCandIDForUps[iMu] = -1;
      L3MuPixCandIDForUps[iMu] = -1;
    }
	
  //0.) check how many L3 muons there are:
  int nMuons = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuCandIDForOnia[iMu] > -1)
      nMuons++;
	
  //   Int_t countCut = 0, countOniaCut = 0;
  //1.) loop over the Pixel tracks
  for (int iP = 0; iP < NohOniaPixel; iP++)
    {
		
      //     countCut = 0;
      //     hEta[histIndex][0][countCut]->Fill(ohOniaPixelEta[iP]);
      //     hPt[histIndex][0][countCut]->Fill(ohOniaPixelPt[iP]);
      //     hHits[histIndex][0][countCut]->Fill(ohOniaPixelHits[iP]);
      //     hNormChi2[histIndex][0][countCut]->Fill(ohOniaPixelNormChi2[iP]);
      //     hDxy[histIndex][0][countCut]->Fill(ohOniaPixelDr[iP]);
      //     hDz[histIndex][0][countCut]->Fill(ohOniaPixelDz[iP]);
		
      //select those that survive the kinematical and
      //topological selection cuts
      if (fabs(ohOniaPixelEta[iP]) > etaPix)
	continue; //eta cut
      if (ohOniaPixelPt[iP] < ptPix)
	continue; //pT cut
      double momThisPix = ohOniaPixelPt[iP] * cosh(ohOniaPixelEta[iP]);
		
      //     hP[histIndex][0][countCut]->Fill(momThisPix);
      //     countCut++;
		
      if (momThisPix < pPix)
	continue; //momentum cut
      if (ohOniaPixelHits[iP] < NHitsPix)
	continue; //min. nb. of hits
      if (ohOniaPixelNormChi2[iP] > normChi2Pix)
	continue; //chi2 cut
      if (fabs(ohOniaPixelDr[iP]) > DxyPix)
	continue; //Dr cut
      if (fabs(ohOniaPixelDz[iP]) > DzPix)
	continue;
		
      //     hEta[histIndex][0][countCut]->Fill(ohOniaPixelEta[iP]);
      //     hPt[histIndex][0][countCut]->Fill(ohOniaPixelPt[iP]);
      //     hHits[histIndex][0][countCut]->Fill(ohOniaPixelHits[iP]);
      //     hNormChi2[histIndex][0][countCut]->Fill(ohOniaPixelNormChi2[iP]);
      //     hDxy[histIndex][0][countCut]->Fill(ohOniaPixelDr[iP]);
      //     hDz[histIndex][0][countCut]->Fill(ohOniaPixelDz[iP]);
      //     hP[histIndex][0][countCut]->Fill(momThisPix);
      //     countCut++;
		
      pix4Mom.SetPtEtaPhiM(
			   ohOniaPixelPt[iP],
			   ohOniaPixelEta[iP],
			   ohOniaPixelPhi[iP],
			   muMass);
      //2.) loop now over all L3 muons and check if they would give a
      //Ups pair:
      for (int iMu = 0; iMu < nMuons; iMu++)
	{
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuCandIDForOnia[iMu]],
			      ohMuL3Eta[L3MuCandIDForOnia[iMu]],
			      ohMuL3Phi[L3MuCandIDForOnia[iMu]],
			      muMass);
	  onia4Mom = pix4Mom + mu4Mom;
	  double oniaMass = onia4Mom.M();
			
	  //       countOniaCut = 0;
	  //       if(oniaMass < 8.0) continue; //Only Ups
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if (oniaMass < massMinPix[0] || oniaMass> massMaxPix[0]) continue; //mass cut
			
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(checkChargePix)
	    if(ohMuL3Chg[iMu] == ohOniaPixelChg[iP]) continue; //charge cut
			
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(fabs(ohMuL3Dz[iMu] - ohOniaPixelDz[iP])> DzMuonPix) continue;
	  //       hOniaEta[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][0][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][0][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][0][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][0][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
			
	  //store the surviving pixel-muon combinations:
	  if(iNPix < 10)
	    {
	      L3PixelCandIDForUps[iNPix] = iP;
	      L3MuPixCandIDForUps[iNPix] = iMu;
	      iNPix++;
	    }
	  //       printf("mu[%d]-pixel[%d] inv. mass %f\n",
	  //           L3MuCandIDForOnia[iMu], iP, oniaMass);
	}
    }
	
  //   hNCand[histIndex][0]->Fill(iNPix);
	
  if(iNPix!=0)
    {
      for(int inP=0;inP<iNPix;inP++)
	{
			
	  //       hPixCandEta[histIndex]->Fill(ohOniaPixelEta[L3PixelCandIDForUps[inP]]);
	  //       hPixCandPt[histIndex]->Fill(ohOniaPixelPt[L3PixelCandIDForUps[inP]]);
	  //       hPixCandP[histIndex]->Fill(ohOniaPixelPt[L3PixelCandIDForUps[inP]] * cosh(ohOniaPixelEta[L3PixelCandIDForUps[inP]]));
			
	  //       if(iNPix>=2){
	  //         for(int jnP=inP+1;jnP<iNPix;jnP++){
	  //            if(inP!=jnP){
	  //              double dEta = fabs(ohOniaPixelEta[L3PixelCandIDForUps[inP]]-ohOniaPixelEta[L3PixelCandIDForUps[jnP]]);
	  //              double dPhi = fabs(ohOniaPixelPhi[L3PixelCandIDForUps[inP]]-ohOniaPixelPhi[L3PixelCandIDForUps[jnP]]);
	  //              if(dPhi>TMath::Pi()) dPhi = 2.0*TMath::Pi()-dPhi;
	  // //              hPixCanddr[histIndex]->Fill(sqrt(pow(dEta,2)+pow(dPhi,2)));
	  //            }
	  //         }
	  //       }
	}
    }
	
  return iNPix;
}

int OHltTree::OpenHltMuTrackPassed_Ups(
				       double ptTrack,
				       double pTrack,
				       double etaTrack,
				       double DxyTrack,
				       double DzTrack,
				       int NHitsTrack,
				       double normChi2Track,
				       double *massMinTrack,
				       double *massMaxTrack,
				       double DzMuonTrack,
				       bool checkChargeTrack,
				       int histIndex)
{
	
  double pixMatchingDeltaR = 0.01;
  const double muMass = 0.105658367;
  TLorentzVector track4Mom, mu4Mom, onia4Mom;
  int iNTrack = 0;
	
  //0.) check how many pixel-muon combinations there are:
  int nComb = 0;
  for (int iMu = 0; iMu < 10; iMu++)
    if (L3MuPixCandIDForUps[iMu] > -1)
      nComb++;
  //   Int_t countCut = 0, countOniaCut = 0;
  //1.) loop over the Tracker tracks
  for (int iT = 0; iT < NohOniaTrack; iT++)
    {
		
      //select those that survive the kinematical and
      //topological selection cuts
      //     countCut++;
      //     hEta[histIndex][1][countCut]->Fill(ohOniaTrackEta[iT]);
      //     hPt[histIndex][1][countCut]->Fill(ohOniaTrackPt[iT]);
      //     hHits[histIndex][1][countCut]->Fill(ohOniaTrackHits[iT]);
      //     hNormChi2[histIndex][1][countCut]->Fill(ohOniaTrackNormChi2[iT]);
      //     hDxy[histIndex][1][countCut]->Fill(ohOniaTrackDr[iT]);
      //     hDz[histIndex][1][countCut]->Fill(ohOniaTrackDz[iT]);
		
      if (fabs(ohOniaTrackEta[iT]) > etaTrack)
	continue; //eta cut
      if (ohOniaTrackPt[iT] < ptTrack)
	continue; //pT cut
      double momThisTrack = ohOniaTrackPt[iT] * cosh(ohOniaTrackEta[iT]);
      //     printf("track[%d] has eta %f, pT %f and mom %f\n",
      //         iT, ohOniaTrackEta[iT], ohOniaTrackPt[iT], momThisTrack);
		
      //     hP[histIndex][1][countCut]->Fill(momThisTrack);
      //     countCut++;
		
      if (momThisTrack < pTrack)
	continue; //momentum cut
      if (ohOniaTrackHits[iT] < NHitsTrack)
	continue; //min. nb. of hits
      if (ohOniaTrackNormChi2[iT] > normChi2Track)
	continue; //chi2 cut
      if (fabs(ohOniaTrackDr[iT]) > DxyTrack)
	continue; //Dr cut
      if (fabs(ohOniaTrackDz[iT]) > DzTrack)
	continue;
		
      //     hEta[histIndex][1][countCut]->Fill(ohOniaTrackEta[iT]);
      //     hPt[histIndex][1][countCut]->Fill(ohOniaTrackPt[iT]);
      //     hHits[histIndex][1][countCut]->Fill(ohOniaTrackHits[iT]);
      //     hNormChi2[histIndex][1][countCut]->Fill(ohOniaTrackNormChi2[iT]);
      //     hDxy[histIndex][1][countCut]->Fill(ohOniaTrackDr[iT]);
      //     hDz[histIndex][1][countCut]->Fill(ohOniaTrackDz[iT]);
      //     hP[histIndex][1][countCut]->Fill(momThisTrack);
		
      //     printf("track %d surviving kinematical pre-selection\n", iT);
      //2.) loop over the pixels candidates to see whether the track
      //under investigation has a match to the pixel track
      bool trackMatched = false;
      for (int iPix = 0; iPix < nComb; iPix++)
	{
			
	  if (trackMatched)
	    break; //in case the track was already matched
	  if (L3PixelCandIDForUps[iPix] < 0)
	    continue; //in case the pixel has been matched to a previous track
			
	  double deltaEta = ohOniaPixelEta[L3PixelCandIDForUps[iPix]]
	    - ohOniaTrackEta[iT];
	  double deltaPhi = ohOniaPixelPhi[L3PixelCandIDForUps[iPix]]
	    - ohOniaTrackPhi[iT];
	  if (deltaPhi>TMath::Pi())
	    deltaPhi = 2.0*TMath::Pi()-deltaPhi;
	  double deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
			
	  //       printf("delta R = %f\n", deltaR);
	  if (deltaR > pixMatchingDeltaR)
	    continue;
	  //       printf("track[%d] and pixel[%d] are compatible (deltaR %f)\n", iT, L3PixelCandIDForUps[iPix], deltaR);
			
	  trackMatched = true;
	  L3PixelCandIDForUps[iPix] = -1; //deactivate this candidate to not match it to any further track
			
	  track4Mom.SetPtEtaPhiM(
				 ohOniaTrackPt[iT],
				 ohOniaTrackEta[iT],
				 ohOniaTrackPhi[iT],
				 muMass);
	  //check if the matched tracker track combined with the
	  //muon gives again an opposite sign onia:
	  mu4Mom.SetPtEtaPhiM(
			      ohMuL3Pt[L3MuPixCandIDForUps[iPix]],
			      ohMuL3Eta[L3MuPixCandIDForUps[iPix]],
			      ohMuL3Phi[L3MuPixCandIDForUps[iPix]],
			      muMass);
	  onia4Mom = track4Mom + mu4Mom;
			
	  double oniaMass = onia4Mom.M();
	  //       printf("mu[%d]-track[%d] inv. mass %f\n",
	  //           L3MuPixCandIDForUps[iPix], iT, oniaMass);
			
	  //       countOniaCut = 0;
	  //       if(oniaMass < 8.0) continue; //Only Ups
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if (oniaMass < massMinTrack[0] || oniaMass> massMaxTrack[0]) continue; //mass cut
			
	  //       printf("surviving: mu[%d]-track[%d] inv. mass %f\n",
	  //           L3MuPixCandIDForUps[iPix], iT, oniaMass);
			
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(checkChargeTrack)
	    if(ohMuL3Chg[L3MuPixCandIDForUps[iPix]] == ohOniaTrackChg[iT]) continue; //charge cut
			
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       countOniaCut++;
			
	  if(fabs(ohMuL3Dz[L3MuPixCandIDForUps[iPix]] - ohOniaTrackDz[iT])> DzMuonTrack) continue; //deltaZ cut
			
	  //       hOniaEta[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta());
	  //       hOniaRap[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity());
	  //       hOniaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Pt());
	  //       hOniaMass[histIndex][1][countOniaCut]->Fill(onia4Mom.M());
	  //       hOniaP[histIndex][1][countOniaCut]->Fill(sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
	  //       hOniaEtaPt[histIndex][1][countOniaCut]->Fill(onia4Mom.Eta(),onia4Mom.Pt());
	  //       hOniaRapP[histIndex][1][countOniaCut]->Fill(onia4Mom.Rapidity(),sqrt(pow(onia4Mom.Px(),2)+pow(onia4Mom.Py(),2)+pow(onia4Mom.Pz(),2)));
			
	  //store the surviving track-muon combinations:
	  if(iNTrack < 10)
	    iNTrack++;
			
	  break; //don't check further pixels... go to next track
	}
    }
	
  //   if(iNTrack > 0)
  //     printf("found %d final candidates!!!\n", iNTrack);
	
  //   hNCand[histIndex][1]->Fill(iNTrack);
  return iNTrack;
}

bool OHltTree::OpenJetID(int jetindex)
{
  if (jetindex>=NohJetCorCal) return false;
  bool jetID = true ; //jetID is true by default
  if (fabs(ohJetCorCalEta[jetindex])< 2.6) {//jetID might be changed to false only for central jets : jetID is a cut only meant for central jets
    jetID =  (ohJetCorCalEMF[jetindex] > 1.0E-6) && (ohJetCorCalEMF[jetindex] < 999.0) && ohJetCorCalN90hits[jetindex]>=2;
  }
  return jetID;
	
}


int OHltTree::OpenHlt1JetPassed(double pt)
{
  int rc = 0;
	
  // Loop over all oh jets 
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i]>pt)
	{ // Jet pT cut
	  rc++;
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1JetPassed(double pt, double etamax)
{
  int rc = 0;
  //ccla
  // Loop over all oh jets 
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i]>pt && fabs(ohJetCalEta[i])<etamax)
	{ // Jet pT cut
	  rc++;
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1JetPassed(
				double pt,
				double etamax,
				double emfmin,
				double emfmax)
{
  int rc = 0;
  //ccla
  // Loop over all oh jets 
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i]>pt && fabs(ohJetCalEta[i])<etamax
	  && ohJetCalEMF[i] > emfmin && ohJetCalEMF[i] < emfmax)
	{ // Jet pT cut
	  rc++;
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1CorJetPassedNoJetID(double pt)
{
  int rc = 0;
	
  // Loop over all oh corrected jets
  for (int i=0; i<NohJetCorCal; i++)
    {
      if ( ohJetCorCalPt[i]>pt)
	{ // Jet pT cut
	  rc++;
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1CorJetPassed(double pt)
{
  int rc = 0;
	
  // Loop over all oh corrected jets
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (OpenJetID(i) && ohJetCorCalPt[i]>pt)
	{ // Jet pT cut
	  rc++;
	}
    }
	
  return rc;
}


int OHltTree::OpenHlt1CorJetPassed(double pt, double etamax)
{
  int rc = 0;
	
  // Loop over all oh corrected jets
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (OpenJetID(i) && ohJetCorCalPt[i]>pt && fabs(ohJetCorCalEta[i])<etamax)
	{ // Jet pT cut
	  rc++;
	}
    }
  return rc;
}

int OHltTree::OpenHlt1PFJetPassed(double pt, double etamax)
{
  int rc = 0; 
         
  // Loop over all oh corrected jets 
  for (int i=0; i<NohPFJet; i++) 
    { 
      if (pfJetPt[i]>pt && fabs(pfJetEta[i])<etamax) 
        { // Jet pT cut 
          rc++; 
        } 
    } 
  return rc; 

}

int OHltTree::OpenHltDiJetAvePassed(double pt)
{
  int rc = 0;
	
  // Loop over all oh jets, select events where the *average* pT of a pair is above threshold
  //std::cout << "FL: NohJetCal = " << NohJetCal << std::endl;
  for (int i=0; i<NohJetCal; i++)
    {
      bool jetID0 = true ; //jetID is true by default
      if (fabs(ohJetCalEta[i])< 2.6) {//jetID might be changed to false only for central jets : jetID is a cut only meant for central jets
	jetID0 =  (ohJetCalEMF[i] > 1.0E-6) && (ohJetCalEMF[i] < 999.0) && ohJetCalN90hits[i]>=2;
      }
      if (!jetID0) continue;
      for (int j=0; j<NohJetCal && j!=i; j++)
	{
	  bool jetID1 = true ; //jetID is true by default
	  if (fabs(ohJetCalEta[i])< 2.6) {//jetID might be changed to false only for central jets : jetID is a cut only meant for central jets
	    jetID1 =  (ohJetCalEMF[i] > 1.0E-6) && (ohJetCalEMF[i] < 999.0) && ohJetCalN90hits[i]>=2;
	  }
	  if (jetID1 && (ohJetCalPt[i]+ohJetCalPt[j])/2.0 > pt)
	    { // Jet pT cut 
	      //      if((ohJetCalE[i]/cosh(ohJetCalEta[i])+ohJetCalE[j]/cosh(ohJetCalEta[j]))/2.0 > pt) {
	      rc++;
	    }
	}
    }
  return rc;
}


int OHltTree::OpenHltCorDiJetAvePassed(double pt)
{
  int rc = 0;
	
  // Loop over all oh jets, select events where the *average* pT of a pair is above threshold
  //std::cout << "FL: NohJetCal = " << NohJetCal << std::endl;
  for (int i=0; i<NohJetCorCal; i++)
    {
      bool jetID0 = true ; //jetID is true by default
      if (fabs(ohJetCorCalEta[i])< 2.6) {//jetID might be changed to false only for central jets : jetID is a cut only meant for central jets
	jetID0 =  (ohJetCorCalEMF[i] > 1.0E-6) && (ohJetCorCalEMF[i] < 999.0) && ohJetCorCalN90hits[i]>=2;
      }
      if (!jetID0) continue;
      for (int j=0; j<NohJetCal && j!=i; j++)
	{
	  bool jetID1 = true ; //jetID is true by default
	  if (fabs(ohJetCorCalEta[i])< 2.6) {//jetID might be changed to false only for central jets : jetID is a cut only meant for central jets
	    jetID1 =  (ohJetCorCalEMF[i] > 1.0E-6) && (ohJetCorCalEMF[i] < 999.0) && ohJetCorCalN90hits[i]>=2;
	  }
	  if (jetID1 && (ohJetCorCalPt[i]+ohJetCorCalPt[j])/2.0 > pt)
	    { // Jet pT cut 
	      //      if((ohJetCalE[i]/cosh(ohJetCalEta[i])+ohJetCalE[j]/cosh(ohJetCalEta[j]))/2.0 > pt) {
	      rc++;
	    }
	}
    }
  return rc;
}


int OHltTree::OpenHltQuadCorJetPassed(double pt)
{
  int njet = 0;
  int rc = 0;
	
  // Loop over all oh jets
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (ohJetCorCalPt[i] > pt && fabs(ohJetCorCalEta[i]) < 5.0)
	{ // Jet pT cut
	  //std::cout << "FL: fires the jet pt cut" << std::endl;
	  njet++;
	}
    }
	
  if (njet >= 4)
    {
      rc = 1;
    }
	
  return rc;
}

int OHltTree::OpenHltQuadJetPassed(double pt)
{
  int njet = 0;
  int rc = 0;
	
  // Loop over all oh jets
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i] > pt && fabs(ohJetCalEta[i]) < 5.0)
	{ // Jet pT cut
	  njet++;
	}
    }
	
  if (njet >= 4)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltQuadJetPassedPlusTauPFId(
					      double pt,
					      double etaJet,
					      double ptTau)
{
  int njet = 0;
  int rc = 0;
  bool foundPFTau = false;
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (ohJetCorCalPt[i] > pt && fabs(ohJetCorCalEta[i]) < etaJet)
	{ // Jet pT cut 
	  njet++;
	  for (int j=0; j<NohpfTau; j++)
	    {
				
	      if (ohpfTauPt[j] > ptTau && ohpfTauLeadTrackPt[j]>= 5
		  && fabs(ohpfTauEta[j]) <2.5 && ohpfTauTrkIso[j] <1
		  && ohpfTauGammaIso[j] <1)
		{
					
		  float deltaEta = ohpfTauEta[j] - ohJetCorCalEta[i];
		  float deltaPhi = ohpfTauPhi[j] - ohJetCorCalPhi[i];
					
		  if (fabs(deltaPhi)>3.141592654)
		    deltaPhi = 6.283185308-fabs(deltaPhi);
					
		  float deltaR = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
					
		  if (deltaR<0.3)
		    {
		      foundPFTau = true;
		    }
		}
	    }
	}
    }
  if (njet >= 4 && foundPFTau == true)
    rc = 1;
  return rc;
}

int OHltTree::OpenL1QuadJet8(double jetPt, double jetEta)
{
  int rc = 0;
	
  for (int i=0; i<NL1CenJet; i++)
    if (L1CenJetEt[i] >= jetPt && fabs(L1CenJetEta[i])<jetEta)
      rc++;
  for (int i=0; i<NL1ForJet; i++)
    if (L1ForJetEt[i] >= jetPt && fabs(L1ForJetEta[i])<jetEta)
      rc++;
  for (int i=0; i<NL1Tau; i++)
    if (L1TauEt [i] >= jetPt && fabs(L1TauEta[i])<jetEta)
      rc++;
	
  return rc;
}

int OHltTree::OpenHltFwdCorJetPassed(double esum)
{
  int rc = 0;
  double gap = 0.;
	
  // Loop over all oh jets, count the sum of energy deposited in HF 
  for (int i=0; i<NohJetCorCal; i++)
    {
      if (((ohJetCorCalEta[i] > 3.0 && ohJetCorCalEta[i] < 5.0)
	   || (ohJetCorCalEta[i] < -3.0 && ohJetCorCalEta[i] > -5.0)))
	{
	  gap+=ohJetCorCalE[i];
	}
    }
	
  // Backward FWD physics logic - we want to select the events *without* large jet energy in HF 
  if (gap < esum)
    rc = 1;
  else
    rc = 0;
	
  return rc;
}

int OHltTree::OpenHltFwdJetPassed(double esum)
{
  int rc = 0;
  double gap = 0.;
	
  // Loop over all oh jets, count the sum of energy deposited in HF 
  for (int i=0; i<NohJetCal; i++)
    {
      if (((ohJetCalEta[i] > 3.0 && ohJetCalEta[i] < 5.0)
	   || (ohJetCalEta[i] < -3.0 && ohJetCalEta[i] > -5.0)))
	{
	  gap+=ohJetCalE[i];
	}
    }
	
  // Backward FWD physics logic - we want to select the events *without* large jet energy in HF 
  if (gap < esum)
    rc = 1;
  else
    rc = 0;
	
  return rc;
}

int OHltTree::OpenHltHTJetNJPassed(
				   double HTthreshold,
				   double jetthreshold,
				   double etamax,
				   int nj)
{
  int rc = 0, njets=0;
  double sumHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold && fabs(ohJetCalEta[i])<etamax)
	{
	  //sumHT+=ohJetCalPt[i];
	  njets++;
	  sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
	}
    }
	
  if (sumHT >= HTthreshold && njets>=nj)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltMHT(double MHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double mhtx=0., mhty=0.;
  for (int i=0; i<NohJetCorCal; ++i)
    {
      if (OpenJetID(i) && ohJetCorCalPt[i] >= jetthreshold && fabs(ohJetCorCalEta[i]) < etathreshold)
	{
	  mhtx-=ohJetCorCalPt[i]*cos(ohJetCorCalPhi[i]);
	  mhty-=ohJetCorCalPt[i]*sin(ohJetCorCalPhi[i]);
	}
    }
  if (sqrt(mhtx*mhtx+mhty*mhty)>MHTthreshold)
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}


int OHltTree::OpenHltPFMHT(double PFMHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double mhtx=0., mhty=0.;
  for (int i=0; i<NohPFJet; ++i)
    {
      if (pfJetPt[i] >= jetthreshold && fabs(pfJetEta[i]) < etathreshold)
	{
	  mhtx-=pfJetPt[i]*cos(pfJetPhi[i]);
	  mhty-=pfJetPt[i]*sin(pfJetPhi[i]);
	}
    }
  if (sqrt(mhtx*mhtx+mhty*mhty)>PFMHTthreshold)
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}


int OHltTree::OpenHltMHTU(double MHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double mhtx=0., mhty=0.;
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold && fabs(ohJetCalEta[i]) < etathreshold)
	{
	  mhtx-=ohJetCalPt[i]*cos(ohJetCalPhi[i]);
	  mhty-=ohJetCalPt[i]*sin(ohJetCalPhi[i]);
	}
    }
  if (sqrt(mhtx*mhtx+mhty*mhty)>MHTthreshold)
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}

int OHltTree::OpenHltPT12U(double PT12threshold, double jetthreshold)
{
  int rc = 0;
  int njets = 0;
  double pt12tx=0., pt12ty=0.;
  for (int i=0; i<NohJetCal; ++i)
    {
      if ((ohJetCalPt[i] >= jetthreshold) && (fabs(ohJetCalEta[i]) <3))
	{
	  njets++;
	  if (njets<3)
	    {
	      pt12tx-=ohJetCalPt[i]*cos(ohJetCalPhi[i]);
	      pt12ty-=ohJetCalPt[i]*sin(ohJetCalPhi[i]);
	    }
	}
    }
  if ((njets >= 2) && (sqrt(pt12tx*pt12tx+pt12ty*pt12ty)>PT12threshold))
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}

int OHltTree::OpenHltPT12(double PT12threshold, double jetthreshold)
{
  int rc = 0;
  int njets = 0;
  double pt12tx=0., pt12ty=0.;
  for (int i=0; i<NohJetCorCal; ++i)
    {
      if ((ohJetCorCalPt[i] >= jetthreshold)
	  && (fabs(ohJetCorCalEta[i]) <3))
	{
	  njets++;
	  if (njets<3)
	    {
	      pt12tx-=ohJetCorCalPt[i]*cos(ohJetCorCalPhi[i]);
	      pt12ty-=ohJetCorCalPt[i]*sin(ohJetCorCalPhi[i]);
	    }
	}
    }
  if ((njets >= 2) && (sqrt(pt12tx*pt12tx+pt12ty*pt12ty)>PT12threshold))
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}

int OHltTree::OpenHltSumHTPassed(double sumHTthreshold, double jetthreshold)
{
  int rc = 0;
  double sumHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold)
	{
	  //sumHT+=ohJetCorCalPt[i];
			
	  sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
	}
    }
	
  if (sumHT >= sumHTthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltSumCorHTPassed(double sumHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double sumHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCorCal; ++i)
    {
      if (ohJetCorCalPt[i] >= jetthreshold && fabs(ohJetCorCalEta[i]) < etathreshold)
	{
	  //sumHT+=recoJetCorCorCalPt[i];
	  sumHT+=ohJetCorCalPt[i];
	  //sumHT+=(ohJetCorCalE[i]/cosh(ohJetCorCalEta[i]));
	}
    }
	
  if (sumHT >= sumHTthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltSumFJCorHTPassed(double sumHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double sumHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCorL1L2L3Cal; ++i)
    {
      if (ohJetCorL1L2L3CalPt[i] >= jetthreshold && fabs(ohJetCorL1L2L3CalEta[i]) < etathreshold)
	{
	  sumHT+=ohJetCorL1L2L3CalPt[i];
	}
    }
	
  if (sumHT >= sumHTthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltSumPFHTPassed(double sumPFHTthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  double sumPFHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohPFJet; ++i)
    {
      if (pfJetPt[i] >= jetthreshold && fabs(pfJetEta[i]) < etathreshold)
	{
	  sumPFHT+=pfJetPt[i];
	}
    }
	
  if (sumPFHT >= sumPFHTthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltMeffU(double Meffthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  //MHT
  double mhtx=0., mhty=0.;
  //HT
  double sumHT = 0.;  
	
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold && fabs(ohJetCalEta[i]) < etathreshold)
	{
	  mhtx-=ohJetCalPt[i]*cos(ohJetCalPhi[i]);
	  mhty-=ohJetCalPt[i]*sin(ohJetCalPhi[i]);
	  sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
	}
    }
	
  if (sqrt(mhtx*mhtx+mhty*mhty)+sumHT>Meffthreshold)
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}

int OHltTree::OpenHltMeff(double Meffthreshold, double jetthreshold, double etathreshold)
{
  int rc = 0;
  //MHT
  double mhtx=0., mhty=0.;
  //HT
  double sumHT = 0.;   
  for (int i=0; i<NohJetCorCal; ++i)
    {
      if (OpenJetID(i) && ohJetCorCalPt[i] >= jetthreshold && fabs(ohJetCorCalEta[i]) < etathreshold)
	{
	  mhtx-=ohJetCorCalPt[i]*cos(ohJetCorCalPhi[i]);
	  mhty-=ohJetCorCalPt[i]*sin(ohJetCorCalPhi[i]);
	  sumHT+=ohJetCorCalPt[i];
	}
    }
	
  if (sqrt(mhtx*mhtx+mhty*mhty)+sumHT>Meffthreshold)
    rc = 1;
  else
    rc = 0;
  //std::cout << "sqrt(mhtx*mhtx+mhty*mhty) = " << sqrt(mhtx*mhtx+mhty*mhty) << std::endl;
  return rc;
}

int OHltTree::OpenHltSumHTPassed(
				 double sumHTthreshold,
				 double jetthreshold,
				 double etajetthreshold)
{
  int rc = 0;
  double sumHT = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold && fabs(ohJetCalEta[i])
	  < etajetthreshold)
	{
	  //sumHT+=ohJetCorCalPt[i];
			
	  sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
	}
    }
	
  if (sumHT >= sumHTthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltSumHTPassed(
				 double sumHTthreshold,
				 double jetthreshold,
				 double etajetthreshold,
				 int Njetthreshold)
{
  int rc = 0;
  double sumHT = 0.;
  int Njet = 0.;
	
  // Loop over all oh jets, sum up the energy  
  for (int i=0; i<NohJetCal; ++i)
    {
      if (ohJetCalPt[i] >= jetthreshold && fabs(ohJetCalEta[i])
	  < etajetthreshold)
	{
	  //sumHT+=ohJetCorCalPt[i];
	  Njet++;
			
	  sumHT+=(ohJetCalE[i]/cosh(ohJetCalEta[i]));
	}
    }
	
  if (sumHT >= sumHTthreshold && Njet >= Njetthreshold)
    rc = 1;
	
  return rc;
}

int OHltTree::OpenHltHT_AlphaT(double HT,double betaT, double Jet){
  int rc = 0;
  double ht = 0.;
  double mhtx = 0., mhty = 0.;
  int nJets = 0;
  double dht = 0.;
  if(betaT < 10.){betaT /= 10.;}
  if(betaT >10.){betaT /= 100.;}
  // cout << "Beta Value is " << betaT << endl;
  for (int i=0; i<NohJetCorCal; i++){
    double mHT = 0.;
    double aT = 0.;
    if(OpenJetID(i) && fabs(ohJetCorCalEta[i]) < 3.0 && (ohJetCorCalE[i]/cosh(ohJetCorCalEta[i])) > Jet ){ // Make ht and mHT for each jet in loop
      ht += (ohJetCorCalE[i]/cosh(ohJetCorCalEta[i])); // HT
      mhtx-=((ohJetCorCalE[i]/cosh(ohJetCorCalEta[i]))*cos(ohJetCorCalPhi[i]));
      mhty-=((ohJetCorCalE[i]/cosh(ohJetCorCalEta[i]))*sin(ohJetCorCalPhi[i]));
      mHT = sqrt(mhty*mhty + mhtx*mhtx); // Make MHT
      nJets++;
      dht += ( nJets < 2 ? (ohJetCorCalE[i]/cosh(ohJetCorCalEta[i])) : -1.* (ohJetCorCalE[i]/cosh(ohJetCorCalEta[i])) ); // Ripped from rob for DHT
      if ( nJets == 2 || nJets == 3 ) {
	aT = ( ht - fabs(dht) ) / ( 2. * sqrt( ( ht*ht ) - ( mHT*mHT ) ) ); // calc alphaT
      } else if ( nJets > 3 ) {
	aT =  ( ht ) / ( 2. * sqrt( ( ht*ht ) - ( mHT*mHT ) ) ); // Calc Beta T if more jets
      }
      if(aT > betaT && ht > HT){
	rc++; // set RC to passed
	return rc; // return RC if any event passes trigger
      }
    }
  }
  return rc; // if no pass this returns zero
}

int OHltTree::OpenHltFatJetPassed(float jetPt, float DR, float DEta, float DiFatJetMass) {

  // list of good jets                                                                                                                                                        
  vector<TLorentzVector> JETS;
  for(int i=0; i<NohJetCorCal; i++) {
    if (fabs(ohJetCorCalEta[i])>=3 || ohJetCorCalPt[i] < 30.)  continue; // require jets with eta<3 

      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(ohJetCorCalPt[i],
		       ohJetCorCalEta[i],
		       ohJetCorCalPhi[i],
		       ohJetCorCalE[i]);
      JETS.push_back(tmp);
  }
  if(JETS.size()<2.) return 0;

  // look for highest-pT two jets                                                                                                                                             
  int iJet1=-99;
  int iJet2=-99;
  int maxPt = 0.;
  for(int i=0; i<int(JETS.size()); i++) {
    if(JETS[i].Pt() > maxPt) {
      maxPt = JETS[i].Pt();
      iJet1 = i;
    }
  }
  maxPt = 0.;
  for(int i=0; i<int(JETS.size()); i++) {
    if(JETS[i].Pt() > maxPt && i != iJet1) {
      maxPt = JETS[i].Pt();
      iJet2 = i;
    }
  }

  // jet recovery
  TLorentzVector JR1 = JETS[iJet1];
  TLorentzVector JR2 = JETS[iJet2];

  //  if(JR1.Pt()<80. || JR2.Pt()<80.) return 0;
  // delta Eta cut
  if(fabs(JR1.Eta()-JR2.Eta())>DEta) return 0;

  for(int i=0; i < int(JETS.size()); i++) {
    if(i == iJet1 || i == iJet2) continue;
    double dR1 = fabs(JETS[iJet1].DeltaR(JETS[i]));
    double dR2 = fabs(JETS[iJet2].DeltaR(JETS[i]));
    if(dR1 < dR2) { // closest to first jet  
      if(dR1 < DR) {
	JR1 = JR1 + JETS[i];
      }
    } else { // closest to second jet   
      if(dR2 < DR) {
	JR2 = JR2 + JETS[i];
      }
    }
  }

    // mass cut  
  TLorentzVector DiJet = JR1+JR2;
  if(DiJet.M()<DiFatJetMass) return 0;
  return 1;
}
				  



int OHltTree::OpenHlt1PixelTrackPassed(float minpt, float minsep, float miniso)
{
  int rc = 0;
	
  // Loop over all oh pixel tracks, check threshold and separation
  for (int i = 0; i < NohPixelTracksL3; i++)
    {
      for (int i = 0; i < NohPixelTracksL3; i++)
	{
	  if (ohPixelTracksL3Pt[i] > minpt)
	    {
	      float closestdr = 999.;
				
	      // Calculate separation from other tracks above threshold
	      for (int j = 0; j < NohPixelTracksL3 && j != i; j++)
		{
		  if (ohPixelTracksL3Pt[j] > minpt)
		    {
		      float dphi = ohPixelTracksL3Phi[i]-ohPixelTracksL3Phi[j];
		      float deta = ohPixelTracksL3Eta[i]-ohPixelTracksL3Eta[j];
		      float dr = sqrt((deta*deta) + (dphi*dphi));
		      if (dr < closestdr)
			closestdr = dr;
		    }
		}
	      if (closestdr > minsep)
		{
		  // Calculate isolation from *all* other tracks without threshold.
		  if (miniso > 0)
		    {
		      int tracksincone = 0;
		      for (int k = 0; k < NohPixelTracksL3 && k != i; k++)
			{
			  float dphi = ohPixelTracksL3Phi[i]-ohPixelTracksL3Phi[k];
			  float deta = ohPixelTracksL3Eta[i]-ohPixelTracksL3Eta[k];
			  float dr = sqrt((deta*deta) + (dphi*dphi));
			  if (dr < miniso)
			    tracksincone++;
			}
		      if (tracksincone == 0)
			rc++;
		    }
		  else
		    rc++;
		}
	    }
	}
    }
	
  return rc;
}

int OHltTree::OpenHlt1L3MuonPassed(double pt, double eta)
{
  //for BTagMu trigger 
	
  int rcL3 = 0;
  // Loop over all oh L3 muons and apply cuts 
  for (int i=0; i<NohMuL3; i++)
    {
      if (ohMuL3Pt[i] > pt && fabs(ohMuL3Eta[i]) < eta)
	{ // L3 pT and eta cut  
	  rcL3++;
	}
    }
	
  return rcL3;
	
}

int OHltTree::OpenHltMhtOverHTPassed(
				     double HardJetThreshold,
				     double HtJetThreshold,
				     double MhtJetThreshold,
				     double MHTovHT,
				     int NoJets)
{
	
  int rc = 0;
  double newHT = 0.;
  double mhtx = 0., mhty = 0.;
  int nJets = 0;
  for (int i=0; i<NohJetCal; i++)
    {
      if (ohJetCalPt[i]>HtJetThreshold)
	{
	  nJets++;
	}
    }
  if ( (ohJetCalE[0]/cosh(ohJetCalEta[0])) > HardJetThreshold
       && (ohJetCalE[1]/cosh(ohJetCalEta[1])) > HardJetThreshold)
    {
      if (nJets > NoJets)
	{
	  //loop over NohJetCal to calculate a new HT for jets above inputJetPt
	  for (int i=0; i<NohJetCal; i++)
	    {
	      if (fabs(ohJetCalEta[i]) > 3.)
		{
		  continue;
		}
	      if ((ohJetCalE[i]/cosh(ohJetCalEta[i])) > MhtJetThreshold)
		{
		  mhtx-=((ohJetCalE[i]/cosh(ohJetCalEta[i]))
			 *cos(ohJetCalPhi[i]));
		  mhty-=((ohJetCalE[i]/cosh(ohJetCalEta[i]))
			 *sin(ohJetCalPhi[i]));
					
		}
	      if ( (ohJetCalE[i]/cosh(ohJetCalEta[i])) > HtJetThreshold)
		{
		  newHT+=((ohJetCalE[i]/cosh(ohJetCalEta[i])));
		}
	    }
	}
      //end calculation of new HT
      if (nJets > 2)
	{
	  if (newHT > 0. && sqrt(mhtx*mhtx+mhty*mhty)/newHT > MHTovHT)
	    {
	      rc++;
	    }
	}
      if (nJets == 2)
	{
	  if ( (newHT - fabs( (ohJetCalE[0]/cosh(ohJetCalEta[0])
			       - (ohJetCalE[1]/cosh(ohJetCalEta[1])))))/(2*sqrt(newHT*newHT
										- (mhtx*mhtx+mhty*mhty))) > sqrt(1.
														 / (4.*(1.-(MHTovHT*MHTovHT)))))
	    {
	      rc++;
	    }
	}
		
    }
	
  // std::cout << "MHT/HT from Jets is " << sqrt(mhtx*mhtx+mhty*mhty)  /  newHT << " HT is: " << newHT << " MHT is: "<<
  // sqrt(mhtx*mhtx+mhty*mhty) << " The Jet threshold is: " << HtJetThreshold << std::endl;
  return rc;
}

int OHltTree::OpenHltMhtOverHTPassedHTthresh(double HT, double MHTovHT)
{
	
  int rc = 0;
  double newHT = 0.;
  double mhtx = 0., mhty = 0.;
  int nJets = 0;
  for (int i=0; i<NohJetCal; i++)
    {
      if ((ohJetCalE[i]/cosh(ohJetCalEta[i]))>30.)
	{
	  nJets++;
	}
    }
  //loop over NohJetCal to calculate a new HT for jets above inputJetPt
  for (int i=0; i<NohJetCal; i++)
    {
      if ((ohJetCalE[i]/cosh(ohJetCalEta[i])) > 20.)
	{
	  mhtx-=((ohJetCalE[i]/cosh(ohJetCalEta[i]))*cos(ohJetCalPhi[i]));
	  mhty-=((ohJetCalE[i]/cosh(ohJetCalEta[i]))*sin(ohJetCalPhi[i]));
	}
      if ((ohJetCalE[i]/cosh(ohJetCalEta[i])) > 30.)
	{
	  newHT+=((ohJetCalE[i]/cosh(ohJetCalEta[i])));
	}
    }
  // }//end calculation of new HT
  if (nJets > 2)
    {
      if (newHT > 0. && sqrt(mhtx*mhtx+mhty*mhty)/newHT > MHTovHT)
	{
	  rc++;
	}
    }
  if (nJets == 2)
    {
      if ( (newHT - fabs( (ohJetCalE[0]/cosh(ohJetCalEta[0])
			   - (ohJetCalE[1]/cosh(ohJetCalEta[1])))))/(2*sqrt(newHT*newHT
									    - (mhtx*mhtx+mhty*mhty))) > sqrt(1. / (4.*(1.-(MHTovHT*MHTovHT)))))
	{
	  rc++;
	}
    }
	
  return rc;
}

vector<int> OHltTree::VectorOpenHlt1PhotonPassed(
						 float Et,
						 int L1iso,
						 float Tiso,
						 float Eiso,
						 float HisoBR,
						 float HisoEC,
						 float HoverEEB,
						 float HoverEEC,
						 float R9,
						 float ClusShapEB,
						 float ClusShapEC)
{
  vector<int> rc;
  // Loop over all oh photons 
  for (int i=0; i<NohPhot; i++)
    {
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < 2.65)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{
		  if (ohPhotTiso[i] < Tiso + 0.002*ohPhotEt[i])
		    {
		      if (ohPhotEiso[i] < Eiso + 0.012*ohPhotEt[i])
			{
			  if ( (TMath::Abs(ohPhotEta[i]) < 1.479
				&& ohPhotHiso[i] < HisoBR + 0.0050*ohPhotEt[i] 
				&& ohPhotClusShap[i] < ClusShapEB
				&& ohPhotR9[i] < R9)
			       ||
			       (1.479 < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < 2.65 
				&& ohPhotHiso[i] < HisoEC + 0.0050*ohPhotEt[i]
				&& ohPhotClusShap[i] < ClusShapEC))
			    {
			      float EcalEnergy = ohPhotEt[i]/(sin(2*atan(exp(0-ohPhotEta[i]))));
			      if ( (TMath::Abs(ohPhotEta[i]) < 1.479
				    && ohPhotHforHoverE[i]/EcalEnergy < HoverEEB)
				   ||
				   ((1.479 < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < 2.65) 
				    && ohPhotHforHoverE[i]/EcalEnergy < HoverEEC))
				if (ohPhotL1Dupl[i] == false) // remove double-counted L1 SCs   
				  rc.push_back(i);
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}


vector<int> OHltTree::VectorOpenHlt1ElectronSamHarperPassed(
							    float Et,
							    int L1iso,
							    float Tisobarrel,
							    float Tisoendcap,
							    float Tisoratiobarrel,
							    float Tisoratioendcap,
							    float HisooverETbarrel,
							    float HisooverETendcap,
							    float EisooverETbarrel,
							    float EisooverETendcap,
							    float hoverebarrel,
							    float hovereendcap,
							    float clusshapebarrel,
							    float clusshapeendcap,
							    float r9barrel,
							    float r9endcap,
							    float detabarrel,
							    float detaendcap,
							    float dphibarrel,
							    float dphiendcap)
{
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  vector<int> rc;
  // Loop over all oh electrons
  for (int i=0; i<NohEle; i++)
    {
      float ohEleHoverE = ohEleHforHoverE[i]/ohEleE[i];
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohEleEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i])
	  < endcapeta)
	isendcap = 1;
		
      if (ohEleEt[i] > Et)
	{
	  if (TMath::Abs(ohEleEta[i]) < endcapeta)
	    {
	      if (ohEleNewSC[i]<=1)
		{
		  if (ohElePixelSeeds[i]>0)
		    {
		      if (ohEleL1iso[i] >= L1iso)
			{ // L1iso is 0 or 1 
			  if (ohEleL1Dupl[i] == false)
			    { // remove double-counted L1 SCs 
			      if ( (isbarrel && ((ohEleHiso[i]/ohEleEt[i])
						 < HisooverETbarrel))
				   || (isendcap && ((ohEleHiso[i]/ohEleEt[i])
						    < HisooverETendcap)))
				{
				  if ( (isbarrel && ((ohEleEiso[i]/ohEleEt[i])
						     < EisooverETbarrel)) || (isendcap
									      && ((ohEleEiso[i]/ohEleEt[i])
										  < EisooverETendcap)))
				    {
				      if ( ((isbarrel) && (ohEleHoverE < hoverebarrel))
					   || ((isendcap) && (ohEleHoverE
							      < hovereendcap)))
					{
					  if ( (isbarrel && (((ohEleTiso[i] < Tisobarrel
							       && ohEleTiso[i] != -999.) || (Tisobarrel
											     == 999.)))) || (isendcap
													     && (((ohEleTiso[i] < Tisoendcap
														   && ohEleTiso[i] != -999.)
														  || (Tisoendcap == 999.)))))
					    {
					      if (((isbarrel) && (ohEleTiso[i]/ohEleEt[i]
								  < Tisoratiobarrel)) || ((isendcap)
											  && (ohEleTiso[i]/ohEleEt[i]
											      < Tisoratioendcap)))
						{
						  if ( (isbarrel && ohEleClusShap[i]
							< clusshapebarrel) || (isendcap
									       && ohEleClusShap[i]
									       < clusshapeendcap))
						    {
						      if ( (isbarrel && ohEleR9[i]
							    < r9barrel) || (isendcap
									    && ohEleR9[i] < r9endcap))
							{
							  if ( (isbarrel
								&& TMath::Abs(ohEleDeta[i])
								< detabarrel)
							       || (isendcap
								   && TMath::Abs(ohEleDeta[i])
								   < detaendcap))
							    {
							      if ( (isbarrel && ohEleDphi[i]
								    < dphibarrel)
								   || (isendcap
								       && ohEleDphi[i]
								       < dphiendcap))
								{
								  rc.push_back(i);
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}

vector<int> OHltTree::VectorOpenHlt1PhotonPassedR9ID(
						     float Et,
						     float R9ID,
						     int L1iso,
						     float Tiso,
						     float Eiso,
						     float HisoBR,
						     float HisoEC,
						     float HoverEEB,
						     float HoverEEC,
						     float R9,
						     float ClusShapEB,
						     float ClusShapEC)
{
  vector<int> rc;
  float barreleta = 1.479;
  float endcapeta = 2.65;
	
  for (int i=0; i<NohPhot; i++)
    {
      int isbarrel = 0;
      int isendcap = 0;
      if (TMath::Abs(ohPhotEta[i]) < barreleta)
	isbarrel = 1;
      if (barreleta < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < endcapeta)
	isendcap = 1;
		
      if (ohPhotEt[i] > Et)
	{
	  if (TMath::Abs(ohPhotEta[i]) < 2.65)
	    {
	      if (ohPhotL1iso[i] >= L1iso)
		{
		  if (ohPhotTiso[i] < Tiso + 0.001*ohPhotEt[i])
		    {
		      if (ohPhotEiso[i] < Eiso + 0.012*ohPhotEt[i])
			{
			  if (((isbarrel)
			       && ohPhotHiso[i] < HisoBR + 0.005*ohPhotEt[i]
			       && ohEleClusShap[i] < ClusShapEB
			       && ohPhotR9[i] < R9)
			      ||
			      ((isendcap) &&
			       ohPhotHiso[i] < HisoEC + 0.005*ohPhotEt[i]
			       && ohEleClusShap[i] < ClusShapEC))
			    {
			      float EcalEnergy = ohPhotEt[i]/(sin(2*atan(exp(0-ohPhotEta[i]))));
			      float ohPhotHoverE = ohPhotHforHoverE[i]/EcalEnergy;
			      if ( ((isbarrel) && (ohPhotHoverE < HoverEEB))
				   || ((isendcap) && (ohPhotHoverE < HoverEEC)))
				if (ohPhotR9ID[i] > R9ID)
				  if (ohPhotL1Dupl[i] == false) // remove double-counted L1 SCs  
				    rc.push_back(i);
			    }
			}
		    }
		}
	    }
	}
    }
	
  return rc;
}


int OHltTree::OpenL1SetSingleJetBit(const float& thresh)
{
	
  // count number of L1 central, forward and tau jets above threshold
	
  int rc=0;
	
  bool CenJet=false, ForJet=false, TauJet=false;
	
  //size_t size = sizeof(L1CenJetEt)/sizeof(*L1CenJetEt);
  const size_t size = 4;
  // cout << thresh << "\t" << size << endl;
	
  int ncenjet=0;
  for (unsigned int i=0; i<size; ++i)
    {
      if (L1CenJetEt[i] >= thresh)
	++ncenjet;
    }
  CenJet=ncenjet>=1;
	
  int nforjet=0;
  for (unsigned int i=0; i<size; ++i)
    {
      if (L1ForJetEt[i] >= thresh)
	++nforjet;
    }
  ForJet=nforjet>=1;
	
  int ntaujet=0;
  for (unsigned int i=0; i<size; ++i)
    {
      if (L1TauEt[i] >= thresh)
	++ntaujet;
    }
  TauJet=ntaujet>=1;
	
  bool L1SingleJet=(CenJet || ForJet || TauJet );
	
  if (L1SingleJet)
    rc=1;
	
  return (rc );
	
}

int OHltTree::OpenHltCleanedTriJetPassed(
    float Et1,
    float Et2,
    float Et3,
    float AbsEtaMax,
    bool cor,
    const std::string& algo,
    float Deta,
    float Mjj,
    bool etaOpposite,
    bool jetID,
    std::vector<int> ohEleIts)
{
   // fill the jet collection
  int NohMyJet = 0;
  float ohMyJetEta[200];
  float ohMyJetPhi[200];
  float ohMyJetEt[200];
  float ohMyJetE[200];
  float ohMyJetEMF[200];
  float ohMyJetN90[200];

  if ( (cor == false) && (algo == "Calo"))
  {
    NohMyJet = NohJetCal;
    for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
    {
      ohMyJetEta[ohMyJetIt] = ohJetCalEta[ohMyJetIt];
      ohMyJetPhi[ohMyJetIt] = ohJetCalPhi[ohMyJetIt];
      ohMyJetEt[ohMyJetIt] = ohJetCalE[ohMyJetIt] * sin(2.*atan(exp(-1.
          *ohJetCalEta[ohMyJetIt])));
      ohMyJetE[ohMyJetIt] = ohJetCalE[ohMyJetIt];
      ohMyJetEMF[ohMyJetIt] = ohJetCalEMF[ohMyJetIt];
      ohMyJetN90[ohMyJetIt] = ohJetCalN90[ohMyJetIt];
    }
  }

  if ( (cor == true) && (algo == "Calo"))
  {
    NohMyJet = NohJetCorCal;
    for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
    {
      ohMyJetEta[ohMyJetIt] = ohJetCorCalEta[ohMyJetIt];
      ohMyJetPhi[ohMyJetIt] = ohJetCorCalPhi[ohMyJetIt];
//          ohMyJetEt[ohMyJetIt] = ohJetCorCalE[ohMyJetIt] * sin(2.*atan(exp(-1.
//                *ohJetCorCalEta[ohMyJetIt])));
      ohMyJetEt[ohMyJetIt] = ohJetCorCalPt[ohMyJetIt];
      ohMyJetE[ohMyJetIt] = ohJetCorCalE[ohMyJetIt];
      ohMyJetEMF[ohMyJetIt] = ohJetCorCalEMF[ohMyJetIt];
      ohMyJetN90[ohMyJetIt] = ohJetCorCalN90[ohMyJetIt];
    }
  }

  if ( (cor == false) && (algo == "PF"))
  {
    NohMyJet = NohPFJet;
    for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
    {
      ohMyJetEta[ohMyJetIt] = pfJetEta[ohMyJetIt];
      ohMyJetPhi[ohMyJetIt] = pfJetPhi[ohMyJetIt];
      ohMyJetEt[ohMyJetIt] = pfJetPt[ohMyJetIt];
      ohMyJetEMF[ohMyJetIt] = -1.;
      ohMyJetN90[ohMyJetIt] = -1.;
    }
  }

   // clean the jet collection from electrons
  int NohCleanedJet = 0;
  float ohCleanedJetEta[200];
  float ohCleanedJetPhi[200];
  float ohCleanedJetEt[200];
  float ohCleanedJetE[200];
  float ohCleanedJetEMF[200];
  float ohCleanedJetN90[200];

  for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
  {
    bool isMatching = false;
    for (unsigned int ohEleIt = 0; ohEleIt < ohEleIts.size(); ++ohEleIt)
      if (deltaR(ohEleEta[ohEleIts.at(ohEleIt)], ohElePhi[ohEleIts.at(ohEleIt)], ohMyJetEta[ohMyJetIt], ohMyJetPhi[ohMyJetIt]) < 0.3)
        isMatching = true;

    if (isMatching == true)
      continue;

    ohCleanedJetEta[NohCleanedJet] = ohMyJetEta[ohMyJetIt];
    ohCleanedJetPhi[NohCleanedJet] = ohMyJetPhi[ohMyJetIt];
    ohCleanedJetEt[NohCleanedJet] = ohMyJetEt[ohMyJetIt];
    ohCleanedJetE[NohCleanedJet] = ohMyJetE[ohMyJetIt];
    ohCleanedJetEMF[NohCleanedJet] = ohMyJetEMF[ohMyJetIt];
    ohCleanedJetN90[NohCleanedJet] = ohMyJetN90[ohMyJetIt];
    ++NohCleanedJet;
  }

   // do the selection
  int rc = 0;
  if (NohCleanedJet < 2)
    return rc;

   // loop on jets
  for (int i = 0; i < NohCleanedJet; ++i)
  {
    if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[i]) < 2.6)
          && (ohCleanedJetEMF[i] < 0.01))
      continue;
    if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[i]) < 2.6)
          && (ohCleanedJetN90[i] < 2))
      continue;

    for (int j = i+1; j < NohCleanedJet; j++)
    {
      if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[j]) < 2.6)
            && (ohCleanedJetEMF[j] < 0.01))
        continue;
      if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[j]) < 2.6)
            && (ohCleanedJetN90[j] < 2))
        continue;

      PtEtaPhiELorentzVector j1(
          ohCleanedJetEt[i],
      ohCleanedJetEta[i],
      ohCleanedJetPhi[i],
      ohCleanedJetE[i]);
      PtEtaPhiELorentzVector j2(
          ohCleanedJetEt[j],
      ohCleanedJetEta[j],
      ohCleanedJetPhi[j],
      ohCleanedJetE[j]);

      if ( (std::max(ohCleanedJetEt[i], ohCleanedJetEt[j]) > Et1)
            && (std::min(ohCleanedJetEt[i], ohCleanedJetEt[j]) > Et2)
            && (deltaEta(ohCleanedJetEta[i], ohCleanedJetEta[j]) > Deta)
            && ((j1+j2).mass() > Mjj) 
            && ( (etaOpposite == true && ohCleanedJetEta[i]*ohCleanedJetEta[j] < 0.) || (etaOpposite == false) )
            && ( (AbsEtaMax > 0 && (TMath::Abs(ohCleanedJetEta[i]) < AbsEtaMax && TMath::Abs(ohCleanedJetEta[j]) < AbsEtaMax)) || AbsEtaMax < 0 ) ){
        for (int k = 0; k < NohCleanedJet; k++)
          
          if ( k != i 
               && k != j 
               && ohCleanedJetEt[k] > Et3 
               && ( (AbsEtaMax > 0 && (TMath::Abs(ohCleanedJetEta[k]) < AbsEtaMax && TMath::Abs(ohCleanedJetEta[k]) < AbsEtaMax)) || AbsEtaMax < 0 ) )
            ++rc;
      }
    }
  }

  return rc;
}



int OHltTree::OpenHltCleanedDiJetPassed(
      float Et1,
      float Et2,
      float AbsEtaMax,
      bool cor,
      const std::string& algo,
      float Deta,
      bool etaOpposite,
      bool jetID,
      std::vector<int> ohEleIts)
{
   // fill the jet collection
   int NohMyJet = 0;
   float ohMyJetEta[200];
   float ohMyJetPhi[200];
   float ohMyJetPt[200];
   float ohMyJetE[200];
   float ohMyJetEMF[200];
   float ohMyJetN90hits[200];

   if ( (cor == false) && (algo == "Calo"))
   {
      NohMyJet = NohJetCal;
      for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
      {
         ohMyJetEta[ohMyJetIt] = ohJetCalEta[ohMyJetIt];
         ohMyJetPhi[ohMyJetIt] = ohJetCalPhi[ohMyJetIt];
         ohMyJetPt[ohMyJetIt] = ohJetCalE[ohMyJetIt] * sin(2.*atan(exp(-1.
               *ohJetCalEta[ohMyJetIt])));
         ohMyJetE[ohMyJetIt] = ohJetCalE[ohMyJetIt];
         ohMyJetEMF[ohMyJetIt] = ohJetCalEMF[ohMyJetIt];
         ohMyJetN90hits[ohMyJetIt] = ohJetCalN90hits[ohMyJetIt];
      }
   }

   if ( (cor == true) && (algo == "Calo"))
   {
      NohMyJet = NohJetCorCal;
      for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
      {
         ohMyJetEta[ohMyJetIt] = ohJetCorCalEta[ohMyJetIt];
         ohMyJetPhi[ohMyJetIt] = ohJetCorCalPhi[ohMyJetIt];
         ohMyJetPt[ohMyJetIt] = ohJetCorCalPt[ohMyJetIt];
         ohMyJetE[ohMyJetIt] = ohJetCorCalE[ohMyJetIt];
         ohMyJetEMF[ohMyJetIt] = ohJetCorCalEMF[ohMyJetIt];
         ohMyJetN90hits[ohMyJetIt] = ohJetCorCalN90hits[ohMyJetIt];
      }
   }

   if ( (cor == false) && (algo == "PF"))
   {
      NohMyJet = NohPFJet;
      for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
      {
         ohMyJetEta[ohMyJetIt] = pfJetEta[ohMyJetIt];
         ohMyJetPhi[ohMyJetIt] = pfJetPhi[ohMyJetIt];
         ohMyJetPt[ohMyJetIt] = pfJetPt[ohMyJetIt];
         ohMyJetEMF[ohMyJetIt] = -1.;
         ohMyJetN90hits[ohMyJetIt] = -1.;
      }
   }

   // clean the jet collection from electrons
   int rc = 0;

   // loop on each good electron
   for (unsigned int ohEleIt = 0; ohEleIt < ohEleIts.size(); ++ohEleIt) {

    int NohCleanedJet = 0;
    float ohCleanedJetEta[200];
    float ohCleanedJetPhi[200];
    float ohCleanedJetEt[200];
    float ohCleanedJetE[200];
    float ohCleanedJetEMF[200];
    float ohCleanedJetN90hits[200];
  
    for (int ohMyJetIt = 0; ohMyJetIt < NohMyJet; ++ohMyJetIt)
    {
        bool isMatching = false;
        // check if this jet is overlapped with the ohEleIt good electron
        if (deltaR(ohEleEta[ohEleIts.at(ohEleIt)], ohElePhi[ohEleIts.at(ohEleIt)], ohMyJetEta[ohMyJetIt], ohMyJetPhi[ohMyJetIt]) < 0.3)
          isMatching = true;
    
        if (isMatching == true)
          continue;
  
        ohCleanedJetEta[NohCleanedJet] = ohMyJetEta[ohMyJetIt];
        ohCleanedJetPhi[NohCleanedJet] = ohMyJetPhi[ohMyJetIt];
        ohCleanedJetEt[NohCleanedJet] = ohMyJetPt[ohMyJetIt];
        ohCleanedJetE[NohCleanedJet] = ohMyJetE[ohMyJetIt];
        ohCleanedJetEMF[NohCleanedJet] = ohMyJetEMF[ohMyJetIt];
        ohCleanedJetN90hits[NohCleanedJet] = ohMyJetN90hits[ohMyJetIt];
        ++NohCleanedJet;
    }
  
    // do the selection
    if (NohCleanedJet < 2)
      continue; 
     
    // loop on jets
    for (int i = 0; i < NohCleanedJet; ++i)
    {
      if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[i]) < 2.6)
            && (ohCleanedJetEMF[i] < 0.01))
        continue;
      if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[i]) < 2.6)
            && (ohCleanedJetN90hits[i] < 2))
        continue;

      for (int j = i+1; j < NohCleanedJet; j++)
      {
        if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[j]) < 2.6)
              && (ohCleanedJetEMF[j] < 0.01))
            continue;
        if ( (jetID == true) && (TMath::Abs(ohCleanedJetEta[j]) < 2.6)
              && (ohCleanedJetN90hits[j] < 2))
            continue;
        
        PtEtaPhiELorentzVector j1(
              ohCleanedJetEt[i],
              ohCleanedJetEta[i],
              ohCleanedJetPhi[i],
              ohCleanedJetE[i]);
        PtEtaPhiELorentzVector j2(
              ohCleanedJetEt[j],
              ohCleanedJetEta[j],
              ohCleanedJetPhi[j],
              ohCleanedJetE[j]);


        if ( (std::max(ohCleanedJetEt[i], ohCleanedJetEt[j]) > Et1)
              && (std::min(ohCleanedJetEt[i], ohCleanedJetEt[j]) > Et2)
              && (deltaEta(ohCleanedJetEta[i], ohCleanedJetEta[j]) > Deta)
              && ( (etaOpposite == true && ohCleanedJetEta[i]*ohCleanedJetEta[j] < 0.) || (etaOpposite == false) )
              && ( (AbsEtaMax > 0 && (TMath::Abs(ohCleanedJetEta[i]) < AbsEtaMax && TMath::Abs(ohCleanedJetEta[j]) < AbsEtaMax)) || AbsEtaMax < 0 ) ) 
            ++rc;

      }
    }// end loop on cleaned jets
   }// end loop on each good electron

   return rc;
}

int OHltTree::OpenHlt1ElectronVbfEleIDPassed(
      float Et,
      float L1SeedEt,
      bool iso,
      int& EtMaxIt,
      std::vector<int>* it,
      bool WP80 = false)
{
   int rc = 0;
   float EtMax = -9999.;
   if (it != NULL)
      it->clear();

  // Loop over all oh electrons  
  for (int i = 0; i < NohEle; ++i)
  {
    // ET/eta/pixelMatch cuts
    if (ohEleEt[i] < Et)
      continue;
    if (TMath::Abs(ohEleEta[i]) > 2.65)
      continue;
    if (ohElePixelSeeds[i] <= 0)
      continue;
    if (ohEleP[i] < 0.)
      continue;

    // EleID
    // R9
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleR9[i] > 999.))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleR9[i] > 999.))
      continue;
    
    if ( WP80 == false ) {
    
    // sigmaietaieta
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleClusShap[i] > 0.011))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleClusShap[i] > 0.031))
      continue;

    // deta
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (TMath::Abs(ohEleDeta[i])
          > 0.008))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (TMath::Abs(ohEleDeta[i])
          > 0.008))
      continue;

    // dphi
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (TMath::Abs(ohEleDphi[i])
          > 0.070))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (TMath::Abs(ohEleDphi[i])
          > 0.050))
      continue;

    // H/E
    if (ohEleHforHoverE[i]/ohEleE[i] > 0.05)
      continue;

    // isolation

    // tracker iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleTiso[i]
            /ohEleEt[i] > 0.125))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleTiso[i]
            /ohEleEt[i] > 0.075))
      continue;

    // ecal iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleEiso[i]
            /ohEleEt[i] > 0.125))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleEiso[i]
            /ohEleEt[i] > 0.075))
      continue;

    // hcal iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleHiso[i]
            /ohEleEt[i] > 0.125))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleHiso[i]
            /ohEleEt[i] > 0.075))
      continue;
    
    }
    
    if ( WP80 == true ) {
    
    // sigmaietaieta
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleClusShap[i] > 0.010))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleClusShap[i] > 0.030))
      continue;

    // deta
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (TMath::Abs(ohEleDeta[i])
          > 0.004))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (TMath::Abs(ohEleDeta[i])
          > 0.007))
      continue;

    // dphi
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && (TMath::Abs(ohEleDphi[i])
          > 0.060))
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && (TMath::Abs(ohEleDphi[i])
          > 0.030))
      continue;

    // H/E
    if ( (TMath::Abs(ohEleEta[i]) < 1.479) && ohEleHforHoverE[i]/ohEleE[i] > 0.04)
      continue;
    if ( (TMath::Abs(ohEleEta[i]) > 1.479) && ohEleHforHoverE[i]/ohEleE[i] > 0.025)
      continue;
    
    // isolation

    // tracker iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleTiso[i]
            /ohEleEt[i] > 0.09))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleTiso[i]
            /ohEleEt[i] > 0.04))
      continue;

    // ecal iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleEiso[i]
            /ohEleEt[i] > 0.07))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleEiso[i]
            /ohEleEt[i] > 0.05))
      continue;

    // hcal iso
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) < 1.479) && (ohEleHiso[i]
            /ohEleEt[i] > 0.1))
      continue;
    if ( (iso == true) && (TMath::Abs(ohEleEta[i]) > 1.479) && (ohEleHiso[i]
            /ohEleEt[i] > 0.025))
      continue;
    
    }

      // L1 matching - isolated

      bool isL1IsolMatched = false;
      for (int j = 0; j < NL1IsolEm; ++j)
      {
         if (L1IsolEmEt[j] > L1SeedEt)
         {
            float etaLow = 0.;
            float etaHigh = 0.;

            if (fabs(ohEleEta[i]) < 1.479)
            {
               etaLow = L1IsolEmEta[j] - 0.522/2.;
               etaHigh = L1IsolEmEta[j] + 0.522/2.;
            }
            else
            {
               etaLow = L1IsolEmEta[j] - 1.0/2.;
               etaHigh = L1IsolEmEta[j] + 1.0/2.;
            }

            float Dphi = deltaPhi(L1IsolEmPhi[j], ohElePhi[i]);

            if (ohEleEta[i] > etaLow && ohEleEta[i] < etaHigh && Dphi < 1.044
                  /2.)
               isL1IsolMatched = true;
         }
      }

      // L1 matching - non isolated
      bool isL1NIsolMatched = false;
      for (int j = 0; j < NL1NIsolEm; ++j)
      {
         if (L1NIsolEmEt[j] > L1SeedEt)
         {
            float etaLow = 0.;
            float etaHigh = 0.;

            if (fabs(ohEleEta[i]) < 1.479)
            {
               etaLow = L1NIsolEmEta[j] - 0.522/2.;
               etaHigh = L1NIsolEmEta[j] + 0.522/2.;
            }
            else
            {
               etaLow = L1NIsolEmEta[j] - 1.0/2.;
               etaHigh = L1NIsolEmEta[j] + 1.0/2.;
            }

            float Dphi = deltaPhi(L1NIsolEmPhi[j], ohElePhi[i]);

            if (ohEleEta[i] > etaLow && ohEleEta[i] < etaHigh && Dphi < 1.044
                  /2.)
               isL1NIsolMatched = true;
         }
      }

      if ( (L1SeedEt > 0.) && (isL1IsolMatched == false) && (isL1NIsolMatched
            == false))
         continue;

      ++rc;

      if (ohEleEt[i] > EtMax)
      {
         EtMaxIt = i;
         EtMax = ohEleEt[i];
      }

      if (it != NULL)
         it->push_back(i);
   }

   return rc;
}



float OHltTree::deltaPhi(const float& phi1, const float& phi2)
{
  float deltaphi = TMath::Abs(phi1 - phi2);
  if (deltaphi > 6.283185308)
    deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654)
    deltaphi = 6.283185308 - deltaphi;
  return deltaphi;
}

float OHltTree::deltaEta(const float& eta1, const float& eta2)
{
  return TMath::Abs(eta1 - eta2);
}

float OHltTree::deltaR(
		       const float& eta1,
		       const float& phi1,
		       const float& eta2,
		       const float& phi2)
{
  return sqrt(deltaEta(eta1, eta2)*deltaEta(eta1, eta2) + deltaPhi(phi1, phi2)
	      *deltaPhi(phi1, phi2));
}

int OHltTree::OpenHltpfMHT(double pfMHTthreshold)
{
  int rc = 0;
  if(pfMHT >= pfMHTthreshold){
    rc = 1;
  }
  return rc;
}

//NOTE: THE MUON AND ELECTRON DR ARE COMMENTED OUT UNTIL THIS IS BETTER UNDERSTOOD
int OHltTree::OpenHltPFTauPassedNoMuonIDNoEleID(float Et,float L25TrkPt, float L3TrkIso, float L3GammaIso,
                                                float mu_ptl1, float mu_ptl2, float mu_ptl3, float mu_dr, float mu_iso,
						float Et_ele, int L1iso,
						float Tisobarrel, float Tisoendcap,
						float Tisoratiobarrel, float Tisoratioendcap,
						float HisooverETbarrel, float HisooverETendcap,
						float EisooverETbarrel, float EisooverETendcap,
						float hoverebarrel, float hovereendcap,
						float clusshapebarrel, float clusshapeendcap,
						float r9barrel, float r9endcap,
						float detabarrel, float detaendcap,
						float dphibarrel, float dphiendcap)
{

  int rc = 0;
  // Loop over all oh pfTaus
  for (int i=0;i < NohpfTau;i++) {
    if (ohpfTauPt[i] >= Et){
      if(fabs(ohpfTauEta[i])<2.5){
	if (ohpfTauLeadTrackPt[i] >= L25TrkPt){
	  if (ohpfTauTrkIso[i] < L3TrkIso){
	    if (ohpfTauGammaIso[i] < L3GammaIso ){
	      // if (OpenHltTauMuonMatching_wMuonID(ohpfTauEta[i], ohpfTauPhi[i], mu_ptl1, mu_ptl2, mu_ptl3, mu_dr, mu_iso) == 0)
	      {
		if (OpenHltTauPFToCaloMatching(ohpfTauEta[i],ohpfTauPhi[i]) == 1){
		  // if (OpenHltTauEleMatching_wEleID(ohpfTauEta[i], ohpfTauPhi[i], Et_ele,L1iso,
		  // 				   Tisobarrel, Tisoendcap,
		  // 				   Tisoratiobarrel, Tisoratioendcap,
		  // 				   HisooverETbarrel, HisooverETendcap,
		  // 				   EisooverETbarrel, EisooverETendcap,
		  // 				   hoverebarrel, hovereendcap,
		  // 				   clusshapebarrel, clusshapeendcap,
		  // 				   r9barrel, r9endcap,
		  // 				   detabarrel, detaendcap,
		  // 				   dphibarrel, dphiendcap) == 0)
		  {

		    rc++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

  }

  return rc;
}


int OHltTree::OpenHltQuadJetCORPassedPlusTauPFIdNewIso(double pt, double etaJet, double ptTau)
{
  int njet=0;
  int rc=0;
  bool foundPFTau=false;
  float deltaR_L1=1000;
  float deltaR_L1Tau=1000;
	
  for(int i=0;i<NrecoJetCorCal;i++){
    if(recoJetCorCalPt[i]>pt&&fabs(recoJetCorCalEta[i])<etaJet)
      njet++;
  }
	
  for(int j=0;j<NohpfTauTightCone;j++){
		
    for(int k=0;k<NL1CenJet;k++){
			
      float deltaETA=ohpfTauTightConeEta[j]-L1CenJetEta[k];
      float deltaPHI=ohpfTauTightConePhi[j]-L1CenJetPhi[k];
      if(fabs(deltaPHI)>3.141592654)deltaPHI=6.283185308-fabs(deltaPHI);
      float deltaR_L1_J=sqrt(pow(deltaETA,2)+pow(deltaPHI,2));
      if(deltaR_L1_J<deltaR_L1)deltaR_L1=deltaR_L1_J;
    }
    for(int s=0;s<NL1Tau;s++){
			
      float deltaETA_Tau=ohpfTauTightConeEta[j]-L1TauEta[s];
      float deltaPHI_Tau=ohpfTauTightConePhi[j]-L1TauPhi[s];
			
      if(fabs(deltaPHI_Tau)>3.141592654)deltaPHI_Tau=6.283185308-fabs(deltaPHI_Tau);
      float deltaR_L1Tau_T=sqrt(pow(deltaETA_Tau,2)+pow(deltaPHI_Tau,2));
      if(deltaR_L1Tau_T<deltaR_L1Tau)deltaR_L1Tau=deltaR_L1Tau_T;
    }
		
    if(deltaR_L1<0.3||deltaR_L1Tau<0.3)
      {
	if(ohpfTauTightConePt[j]>ptTau  && ohpfTauTightConeLeadTrackPt[j]>=5 && fabs(ohpfTauTightConeEta[j])<2.5 && ohpfTauTightConeTrkIso[j]<1.0 &&
	   ohpfTauTightConeGammaIso[j]<1.5)
	  {
	    foundPFTau=true;
	  }
      }
  }
	
	
  if(njet>=4&&foundPFTau==true) rc=1;
	
  return rc;
}

int OHltTree::OpenHlt2MuonOSMassPassed(
				       double ptl1,
				       double ptl2,
				       double ptl3,
				       double dr,
				       int iso,
				       double masslow,
				       double masshigh)
{
  int rc = 0;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TLorentzVector diMu;
  const double muMass = 0.105658367;
  
  for (int i=0; i<NohMuL3; i++)
    {
      for (int j=i+1; j<NohMuL3; j++)
	{
	  if (ohMuL3Pt[i] > ptl3 && ohMuL3Pt[j] > ptl3)
	    { // L3 pT cut
	      if (ohMuL3Dr[i] < dr && ohMuL3Dr[j] < dr)
		{ // L3 DR cut
		  if ((ohMuL3Chg[i] * ohMuL3Chg[j]) < 0)
		    { // opposite charge
		      
		      mu1.SetPtEtaPhiM(
				       ohMuL3Pt[i],
				       ohMuL3Eta[i],
				       ohMuL3Phi[i],
				       muMass);
		      mu2.SetPtEtaPhiM(
				       ohMuL3Pt[j],
				       ohMuL3Eta[j],
				       ohMuL3Phi[j],
				       muMass);
		      diMu = mu1 + mu2;
		      float diMuMass = diMu.M();
		      
		      if ((diMuMass <= masshigh) && (diMuMass >= masslow))
			{
			  int l2match1 = ohMuL3L2idx[i];
			  int l2match2 = ohMuL3L2idx[j];
			  
			  if ( (ohMuL2Pt[l2match1] > ptl2)
			       && (ohMuL2Pt[l2match2] > ptl2))
			    { // L2 pT cut
			      rc++;
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}

int OHltTree::OpenHlt2MuonOSMassDCAPassed(
					  double ptl1,
					  double ptl2,
					  double ptl3,
					  double dr,
					  int iso,
					  double masslow,
					  double masshigh,
					  double etal3,
					  double dimupt,
					  double dca)
{
  int rc = 0;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TLorentzVector diMu;
  const double muMass = 0.105658367;
  
  for (int i=0; i<NohDiMu; i++) 
    {
      if (ohDiMuDCA[i]<dca)
	{
	  if (ohMuL3Pt[ohDiMu1st[i]] > ptl3 && ohMuL3Pt[ohDiMu2nd[i]] > ptl3)
	    { // L3 pT cut
	      if (fabs(ohMuL3Eta[ohDiMu1st[i]]) < etal3 && fabs(ohMuL3Eta[ohDiMu2nd[i]]) < etal3)
		{ // L3 eta cut
		  if (ohMuL3Dr[ohDiMu1st[i]] < dr && ohMuL3Dr[ohDiMu2nd[i]] < dr)
		    { // L3 DR cut
		      if ((ohMuL3Chg[ohDiMu1st[i]] * ohMuL3Chg[ohDiMu2nd[i]]) < 0)
			{ // opposite charge
			  
			  mu1.SetPtEtaPhiM(
					   ohMuL3Pt[ohDiMu1st[i]],
					   ohMuL3Eta[ohDiMu1st[i]],
					   ohMuL3Phi[ohDiMu1st[i]],
					   muMass);
			  mu2.SetPtEtaPhiM(
					   ohMuL3Pt[ohDiMu2nd[i]],
					   ohMuL3Eta[ohDiMu2nd[i]],
					   ohMuL3Phi[ohDiMu2nd[i]],
					   muMass);
			  diMu = mu1 + mu2;
			  float diMuMass = diMu.M();
			  float diMuPt = diMu.Pt();
			  
			  if ((diMuMass <= masshigh) && (diMuMass >= masslow) && diMuPt > dimupt)
			    {
			      int l2match1 = ohMuL3L2idx[ohDiMu1st[i]];
			      int l2match2 = ohMuL3L2idx[ohDiMu2nd[i]];
			      
			      if ( (ohMuL2Pt[l2match1] > ptl2)
				   && (ohMuL2Pt[l2match2] > ptl2))
				{ // L2 pT cut
				  rc++;
				}
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}

int OHltTree::OpenHlt2MuonOSMassVtxPassed(
					  double ptl1,
					  double ptl2,
					  double ptl3,
					  double dr,
					  int iso,
					  double masslow,
					  double masshigh,
					  double etal3,
					  double dimupt,
					  double dca,
					  double chi2,
					  double cos,
					  double lxysig)
{
  int rc = 0;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TLorentzVector diMu;
  const double muMass = 0.105658367;
  
  for (int i=0; i<NohDiMu; i++) 
    {
      if (ohDiMuDCA[i]<dca)
	{
	  for (int j=0 ; j<NohDiMuVtx; j++)
	    {
	      if ((ohDiMu1st[i]==ohDiMuVtx1st[j] && ohDiMu2nd[i]==ohDiMuVtx2nd[j]) || (ohDiMu1st[i]==ohDiMuVtx2nd[j] && ohDiMu2nd[i]==ohDiMuVtx1st[j]))
		{
		  if (ohDiMuVtxROverSig[j]>=lxysig && ohDiMuVtxCosAlpha[j]>cos && ohDiMuVtxChi2[j]<chi2)
		    {
		      if (ohMuL3Pt[ohDiMu1st[i]] > ptl3 && ohMuL3Pt[ohDiMu2nd[i]] > ptl3)
			{ // L3 pT cut
			  if (fabs(ohMuL3Eta[ohDiMu1st[i]]) < etal3 && fabs(ohMuL3Eta[ohDiMu2nd[i]]) < etal3)
			    { // L3 eta cut
			      if (ohMuL3Dr[ohDiMu1st[i]] < dr && ohMuL3Dr[ohDiMu2nd[i]] < dr)
				{ // L3 DR cut
				  if ((ohMuL3Chg[ohDiMu1st[i]] * ohMuL3Chg[ohDiMu2nd[i]]) < 0)
				    { // opposite charge
				      mu1.SetPtEtaPhiM(
						       ohMuL3Pt[ohDiMu1st[i]],
						       ohMuL3Eta[ohDiMu1st[i]],
						       ohMuL3Phi[ohDiMu1st[i]],
						       muMass);
				      mu2.SetPtEtaPhiM(
						       ohMuL3Pt[ohDiMu2nd[i]],
						       ohMuL3Eta[ohDiMu2nd[i]],
						       ohMuL3Phi[ohDiMu2nd[i]],
						       muMass);
				      diMu = mu1 + mu2;
				      float diMuMass = diMu.M();
				      float diMuPt = diMu.Pt();
				      

				      if ((diMuMass <= masshigh) && (diMuMass >= masslow) && diMuPt > dimupt)
					{
					  
					  int l2match1 = ohMuL3L2idx[ohDiMu1st[i]];
					  int l2match2 = ohMuL3L2idx[ohDiMu2nd[i]];
					  
					  if ( (ohMuL2Pt[l2match1] > ptl2)
					       && (ohMuL2Pt[l2match2] > ptl2))
					    { // L2 pT cut
					      rc++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}

int OHltTree::OpenHlt2MuonOSMassVtxHiggsPassed(
					  double ptl1,
					  double ptl2,
					  double ptl3,
					  double dr,
					  int iso,
					  double masslow,
					  double masshigh,
					  double etal3,
					  double dimupt,
					  double dca,
					  double chi2,
					  double cos,
					  double lxysig)
{
  int rc = 0;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TLorentzVector diMu;
  const double muMass = 0.105658367;
  
  for (int i=0; i<NohDiMu; i++) 
    {
      if (ohDiMuDCA[i]<dca)
	{
	  for (int j=0 ; j<NohDiMuVtx; j++)
	    {
	      if ((ohDiMu1st[i]==ohDiMuVtx1st[j] && ohDiMu2nd[i]==ohDiMuVtx2nd[j]) || (ohDiMu1st[i]==ohDiMuVtx2nd[j] && ohDiMu2nd[i]==ohDiMuVtx1st[j]))
		{
		  if (ohDiMuVtxROverSig[j]<=lxysig && ohDiMuVtxCosAlpha[j]>cos && ohDiMuVtxChi2[j]<chi2)
		    {
		      if (ohMuL3Pt[ohDiMu1st[i]] > ptl3 && ohMuL3Pt[ohDiMu2nd[i]] > ptl3)
			{ // L3 pT cut
			  if (fabs(ohMuL3Eta[ohDiMu1st[i]]) < etal3 && fabs(ohMuL3Eta[ohDiMu2nd[i]]) < etal3)
			    { // L3 eta cut
			      if (ohMuL3Dr[ohDiMu1st[i]] < dr && ohMuL3Dr[ohDiMu2nd[i]] < dr)
				{ // L3 DR cut
				  if (ohMuL3Iso[i] >= iso)
				    {//L3 iso cut
				      if ((ohMuL3Chg[ohDiMu1st[i]] * ohMuL3Chg[ohDiMu2nd[i]]) < 0)
					{ // opposite charge
					  mu1.SetPtEtaPhiM(
							   ohMuL3Pt[ohDiMu1st[i]],
							   ohMuL3Eta[ohDiMu1st[i]],
							   ohMuL3Phi[ohDiMu1st[i]],
							   muMass);
					  mu2.SetPtEtaPhiM(
							   ohMuL3Pt[ohDiMu2nd[i]],
							   ohMuL3Eta[ohDiMu2nd[i]],
							   ohMuL3Phi[ohDiMu2nd[i]],
							   muMass);
					  diMu = mu1 + mu2;
					  float diMuMass = diMu.M();
					  float diMuPt = diMu.Pt();
					  
					  
					  if ((diMuMass <= masshigh) && (diMuMass >= masslow) && diMuPt > dimupt)
					    {
					      
					      int l2match1 = ohMuL3L2idx[ohDiMu1st[i]];
					      int l2match2 = ohMuL3L2idx[ohDiMu2nd[i]];
					      
					      if ( (ohMuL2Pt[l2match1] > ptl2)
						   && (ohMuL2Pt[l2match2] > ptl2))
						{ // L2 pT cut
						  rc++;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}

int OHltTree::OpenHlt3MuonOSMassVtxPassed(
					  double ptl1,
					  double ptl2,
					  double ptl3,
					  double dr,
					  int iso,
					  double masslow,
					  double masshigh,
					  double etal3,
					  double dimupt,
					  double dca,
					  double chi2,
					  double cos,
					  double lxysig)
{
  int rc = 0;
  TLorentzVector mu1;
  TLorentzVector mu2;
  TLorentzVector diMu;
  const double muMass = 0.105658367;
  
  for (int i=0; i<NohDiMu; i++) 
    {
      if (ohDiMuDCA[i]<dca)
	{
	  for (int j=0 ; j<NohDiMuVtx; j++)
	    {
	      if ((ohDiMu1st[i]==ohDiMuVtx1st[j] && ohDiMu2nd[i]==ohDiMuVtx2nd[j]) || (ohDiMu1st[i]==ohDiMuVtx2nd[j] && ohDiMu2nd[i]==ohDiMuVtx1st[j]))
		{
		  if (ohDiMuVtxROverSig[j]>lxysig && ohDiMuVtxCosAlpha[j]>cos && ohDiMuVtxChi2[j]<chi2)
		    {
		      if (ohMuL3Pt[ohDiMu1st[i]] > ptl3 && ohMuL3Pt[ohDiMu2nd[i]] > ptl3)
			{ // L3 pT cut
			  if (fabs(ohMuL3Eta[ohDiMu1st[i]]) < etal3 && fabs(ohMuL3Eta[ohDiMu2nd[i]]) < etal3)
			    { // L3 eta cut
			      if (ohMuL3Dr[ohDiMu1st[i]] < dr && ohMuL3Dr[ohDiMu2nd[i]] < dr)
				{ // L3 DR cut
				  if ((ohMuL3Chg[ohDiMu1st[i]] * ohMuL3Chg[ohDiMu2nd[i]]) < 0)
				    { // opposite charge
				      
				      mu1.SetPtEtaPhiM(
						       ohMuL3Pt[ohDiMu1st[i]],
						       ohMuL3Eta[ohDiMu1st[i]],
						       ohMuL3Phi[ohDiMu1st[i]],
						       muMass);
				      mu2.SetPtEtaPhiM(
						       ohMuL3Pt[ohDiMu2nd[i]],
						       ohMuL3Eta[ohDiMu2nd[i]],
						       ohMuL3Phi[ohDiMu2nd[i]],
						       muMass);
				      diMu = mu1 + mu2;
				      float diMuMass = diMu.M();
				      float diMuPt = diMu.Pt();
				      
				      if ((diMuMass <= masshigh) && (diMuMass >= masslow) && diMuPt > dimupt)
					{
					  for (int k=0 ; k<NohMuL3 ; k++)
					    {
					      if (k!=ohDiMu1st[i] && k!=ohDiMu2nd[i] && ohMuL3Pt[k] > ptl3 && fabs(ohMuL3Eta[k]) < etal3)
  						{
						  int l2match1 = ohMuL3L2idx[ohDiMu1st[i]];
						  int l2match2 = ohMuL3L2idx[ohDiMu2nd[i]];
						  int l2match3 = ohMuL3L2idx[k];
						  
						  if ( (ohMuL2Pt[l2match1] > ptl2)
						       && (ohMuL2Pt[l2match2] > ptl2)
						       && (ohMuL2Pt[l2match3] > ptl2))
						    { // L2 pT cut
						      rc++;
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  return rc;
}


bool OHltTree::OpenHltNJetPtPassed(int N, const double& pt)
{
    Int_t NpassPt= 0;
    for (int i= 0; i < NohJetCorCal; ++i)
    {
        if (ohJetCorCalPt[i] >= pt) NpassPt++;
    }
    return NpassPt >= N;
}

bool OHltTree::OpenHltNPFJetPassed(const int N, const double& pt, const double& eta)
{
    Int_t Npass= 0;
    for (int i= 0; i < NohPFJet; ++i)
    {
      if (pfJetPt[i] >= pt && abs(pfJetEta[i]) < eta) 
	Npass++;
      if (Npass >= N)
	break;
    }
    return Npass >= N;
}

bool OHltTree::OpenHltNTowCalEtPassed(int N, const double& Et)
{
    Int_t NpassEt= 0;
    for (int i= 0; i < NrecoTowCal; ++i)
    {
        if (recoTowEt[i] >= Et) NpassEt++;
    }
    return NpassEt >= N;
}


bool OHltTree::OpenHltInvMassCutMu(int nMu, const float& invMassCut){

  float invMassMin = 0.;

  for (int i=0; i<nMu; i++)
    {
      TLorentzVector mu_i;
      mu_i.SetPtEtaPhiM(ohMuL3Pt[i], ohMuL3Eta[i], ohMuL3Phi[i], 0.105);

      for (int j = i+1 ; j <nMu ; j++){

	TLorentzVector mu_j;
	mu_j.SetPtEtaPhiM(ohMuL3Pt[j], ohMuL3Eta[j], ohMuL3Phi[j], 0.105);

	float invMass = (mu_i + mu_j).M();
	if (invMass > invMassMin){
	  invMassMin = invMass;
	}	  
      }
    }

  return (invMassMin > invMassCut );

}

bool OHltTree::OpenHltInvMassCutEle(int nEle, const float& invMassCut){

  float invMassMin = 0.;

  for (int i=0; i<nEle; i++)
    {
      TLorentzVector ele_i;
      ele_i.SetPtEtaPhiM(ohEleEt[i], ohEleEta[i], ohElePhi[i], 0.0005);

      for (int j = i+1 ; j <nEle ; j++){

	TLorentzVector ele_j;
	ele_j.SetPtEtaPhiM(ohEleEt[j], ohEleEta[j], ohElePhi[j], 0.0005);

	float invMass = (ele_i + ele_j).M();
	if (invMass > invMassMin){
	  invMassMin = invMass;
	}	  
      }
    }

  return (invMassMin > invMassCut );

}

bool OHltTree::OpenHltInvMassCutEleMu(int nEle, int nMu, const float& invMassCut){

  float invMassMin = 0.;

  for (int i=0; i<nEle; i++)
    {
      TLorentzVector ele_i;
      ele_i.SetPtEtaPhiM(ohEleEt[i], ohEleEta[i], ohElePhi[i], 0.0005);

      for (int j = 0 ; j <nMu ; j++){

	TLorentzVector mu_j;
	mu_j.SetPtEtaPhiM(ohMuL3Pt[j], ohMuL3Eta[j], ohMuL3Phi[j], 0.105);

	float invMass = (ele_i + mu_j).M();
	if (invMass > invMassMin){
	  invMassMin = invMass;
	}	  
      }
    }

  return (invMassMin > invMassCut );

}
