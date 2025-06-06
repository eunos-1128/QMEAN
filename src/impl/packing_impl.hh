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

#ifndef PACKING_IMPL_HH
#define PACKING_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <qmean/atom_types.hh>
#include <limits>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN PackingOpts {

  PackingOpts () { }

  PackingOpts(Real co, int mc, int bs): cutoff(co), max_counts(mc), bin_size(bs), sigma(0.02){ 
    
    if(cutoff <= 0.0){
      throw std::runtime_error("Cutoff must be larger 0.0!");
    }

    if(max_counts <= 0){
      throw std::runtime_error("Max counts must be larger 0!");
    }

    if(bin_size != 1){
      throw std::runtime_error("Bin size must be one!");
    }

    if(sigma <= 0.0){
      throw std::runtime_error("Sigma parameter must be larger 0.0!"); 
    }
  }
  PackingOpts(Real co, int mc, int bs, Real s): cutoff(co), max_counts(mc), bin_size(bs), sigma(s){ 

    if(cutoff <= 0.0){
      throw std::runtime_error("Cutoff must be larger 0.0!");
    }

    if(max_counts <= 0){
      throw std::runtime_error("Max counts , must be larger 0!");
    }

    if(bin_size != 1){
      throw std::runtime_error("Bin size must be one!");
    }

    if(sigma <= 0.0){
      throw std::runtime_error("Sigma parameter must be larger 0.0!"); 
    }
  }

  inline bool operator == (const PackingOpts& opts) const {
    return max_counts == opts.max_counts && cutoff == opts.cutoff &&
           bin_size == opts.bin_size && sigma == opts.sigma;
  }

  inline bool operator != (const PackingOpts& opts) const {
    return !(*this == opts);
  }


  Real cutoff;
  int max_counts;
  int bin_size;
  Real sigma;
  template <typename DS>
  void Serialize(DS& ds){
    ds & cutoff;
    ds & max_counts;
    ds & bin_size;
    ds & sigma;
  }
};

class DLLEXPORT_QMEAN PackingPotentialImpl : public ost::mol::EntityVisitor {
public:

  PackingPotentialImpl() { }

  PackingPotentialImpl(PackingOpts& opts):
    opts_(opts) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& bins)=0;

protected:

  int GetBin(int count);

  PackingOpts opts_;
  ost::mol::EntityView  env_;
};

}} //namespace




#endif // PACKING_IMPL_HH

