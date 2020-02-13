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

#ifndef INTERACTION_IMPL_HH
#define INTERACTION_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <qmean/atom_types.hh>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN InteractionOpts{
InteractionOpts() { }

InteractionOpts(Real l_cutoff, Real u_cutoff, int nob, int ssep, Real s=0.02):
             lower_cutoff(l_cutoff), upper_cutoff(u_cutoff), number_of_bins(nob),
             sequence_sep(ssep), sigma(s){ 

  if(lower_cutoff >= upper_cutoff){
    throw std::runtime_error("Lower cutoff must be smaller than upper cutoff!");
  }
  
  if(number_of_bins <= 0){
    throw std::runtime_error("Number of distance bins must be larger than 0!");
  }

  if(sequence_sep < 0){
    throw std::runtime_error("Sequence separation must be larger or equal 0!");
  }

  if(sigma <= 0.0){
    throw std::runtime_error("Sigma parameter must be larger 0.0!");
  }
}

public:
  Real lower_cutoff;
  Real upper_cutoff;
  int number_of_bins;
  int sequence_sep;
  Real sigma;

  inline bool operator == (const InteractionOpts& opts) const {
    return lower_cutoff == opts.lower_cutoff && upper_cutoff == opts.upper_cutoff &&
           number_of_bins == opts.number_of_bins && sequence_sep == opts.sequence_sep &&
           sigma == opts.sigma;
  }

  inline bool operator != (const InteractionOpts& opts) const {
    return !(*this == opts);
  }

  template <typename DS>
  void Serialize(DS& ds)
  {
    ds & lower_cutoff;
    ds & upper_cutoff;
    ds & number_of_bins;
    ds & sequence_sep;
    ds & sigma;
  }
};

class DLLEXPORT_QMEAN InteractionPotentialImpl : public ost::mol::EntityVisitor {
public:

  InteractionPotentialImpl() { }

  InteractionPotentialImpl(InteractionOpts& opts):
    opts_(opts) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual bool VisitAtom(const ost::mol::AtomHandle& at);

  virtual void OnInteraction(atom::ChemType a,
                             atom::ChemType b,
                             Real dist)=0;

protected:
  InteractionOpts opts_;
  ost::mol::EntityView  env_;

private:
  int actual_residue_number_;
  ost::conop::AminoAcid actual_aa_;
  String actual_chain_name_;

};

}} //namespace




#endif // MREDUCED_IMPL_HH

