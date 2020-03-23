// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
// Biozentrum - University of Basel
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

//What you get at the end is the logodd of the probability of the tuple (dssp_state, psipred_state, psipred_conf)
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

  Real GetOverallP(char d, char psipred, int psiconf) const {
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

