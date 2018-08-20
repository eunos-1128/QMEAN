// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

