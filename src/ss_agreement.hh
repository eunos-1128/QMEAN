#ifndef QMEAN_SSAGREEMENT_HH
#define QMEAN_SSAGREEMENT_HH


#include <ost/base.hh>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

namespace qmean {

class SSAgreement;
typedef boost::shared_ptr<SSAgreement> SSAgreementPtr;


int PsiPredToIndex(char pp);
int DSSPToIndex(char dssp);

//Holds dssp, psipred and confidence value distributions, such that the likelihood can be 
//calculated to observe a certain dssp state given a psipred prediction with its confidence
//value. The data has been extracted from a representative list of CATH groups with a seqid 
//threshold of 0.35. Structures with resolution below 2 A have been filtered out.
//
//Wat you get at the end is the logodd of the probability of the tuple (dssp_state, psipred_state, psipred_conf)
//and the product of the probabilities for this dssp state and the pair of psipred_values
// => log(P(dssp,psipred,conf)/(P(dssp)*P(psipred))

class SSAgreement{
public:

  SSAgreement();

  Real GetPsiPredP(char pred, int conf) const {
    return psipredP_[conf][PsiPredToIndex(pred)];
  }

  Real GetDSSPP(char d) const {
    return dsspP_[DSSPToIndex(d)];
  }

  Real GetConditionalP(char d, char psipred, int psiconf) const {
    if(psiconf<0 || psiconf>9){
      throw std::invalid_argument("confidence value must be in [0,9]");
    }
    return overallP_[DSSPToIndex(d)][PsiPredToIndex(psipred)][psiconf];
  }

  Real LogOddScore(char d, char psipred, int psiconf) const {
    if(psiconf<0 || psiconf>9){
      throw std::invalid_argument("confidence value must be in [0,9]");
    }
    int idx_dssp = DSSPToIndex(d);
    int idx_psipred = PsiPredToIndex(psipred);
    return log(overallP_[idx_dssp][idx_psipred][psiconf]/(psipredP_[psiconf][idx_psipred]*dsspP_[idx_dssp]));
  }

private:
  Real psipredP_[10][3];
  Real dsspP_[7];
  Real overallP_[7][3][10];
};

}

#endif

