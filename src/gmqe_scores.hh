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

#ifndef GMQE_SCORES_HH
#define GMQE_SCORES_HH

#include <qmean/disco.hh>
#include <qmean/potential_base.hh>
#include <qmean/cb_packing_potential.hh>
#include <qmean/reduced_potential.hh>
#include <qmean/torsion_potential.hh>
#include <qmean/cbeta_potential.hh>
#include <qmean/ss_agreement.hh>


namespace qmean {



class GMQEScoreCalculator;
typedef boost::shared_ptr<GMQEScoreCalculator> GMQEScoreCalculatorPtr;

class GMQEScoreCalculator {

public:

  GMQEScoreCalculator(PotentialContainerPtr potentials,
                      DisCoContainerPtr disco_container,
                      const ost::seq::SequenceHandle& seqres,
                      const ost::seq::SequenceHandle& psipred_pred,
                      const std::vector<int>& psipred_cfi);


  void Eval(const geom::Vec3List& n_positions,
            const geom::Vec3List& ca_positions,
            const geom::Vec3List& c_positions,
            const geom::Vec3List& cb_positions,
            const ost::seq::SequenceHandle& dssp_states,
            const std::vector<int>& residue_numbers,
            Real& dist_const, Real& reduced, 
            Real& cb_packing, Real& torsion,
            Real& ss_agreement, Real& coverage) const;

private:

  std::vector<CBPackingEnergies*> cb_packing_energies_;
  std::vector<ReducedEnergies*> reduced_energies_;
  std::vector<TorsionEnergies*> torsion_energies_;
  DisCoContainerPtr disco_container_;
  std::vector<ost::conop::AminoAcid> amino_acids_;
  std::vector<char> psipred_pred_;
  std::vector<int> psipred_cfi_;
  Real reduced_cutoff_;
  Real reduced_squared_cutoff_;
  int reduced_seq_sep_;
  Real cb_packing_squared_cutoff_;  
  int cb_packing_max_count_;
  Real disco_cutoff_;
  Real disco_squared_cutoff_;
  Real disco_bin_size_;
  PotentialContainerPtr potential_container_;
  SSAgreement ss_agreement_;
};

} // ns

#endif
