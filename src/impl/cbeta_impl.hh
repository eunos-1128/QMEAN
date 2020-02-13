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

#ifndef CBETA_IMPL_HH
#define CBETA_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <ost/mol/spatial_organizer.hh>
#include <ost/mol/alg/construct_cbeta.hh>

namespace qmean { namespace impl {

  
struct DLLEXPORT CBetaOpts{
CBetaOpts() { }

CBetaOpts(Real l_cutoff, Real u_cutoff, int nob, int ssep, Real s=0.02):
          lower_cutoff(l_cutoff), upper_cutoff(u_cutoff), number_of_bins(nob),
          sequence_sep(ssep), sigma(s)
  { 
    if(lower_cutoff >= upper_cutoff){
      throw std::runtime_error("Lower cutoff must be smaller than upper cutoff!");
    }

    if(lower_cutoff < 0.0){
      throw std::runtime_error("Lower cutoff must be greater or equal 0.0!");
    }

    if(sigma <= 0.0){
      throw std::runtime_error("Sigma must be larger larger than 0.0!");
    }
    if(sequence_sep < 0){
      throw std::runtime_error("Sequence separation must be larger or equal 0!");
    }
  }

public:
  Real lower_cutoff;
  Real upper_cutoff;
  int number_of_bins;
  int sequence_sep;
  Real sigma;

  inline bool operator == (const CBetaOpts& opts) const {
    return lower_cutoff == opts.lower_cutoff && upper_cutoff == opts.upper_cutoff &&
           number_of_bins == opts.number_of_bins && sequence_sep == opts.sequence_sep &&
           sigma == opts.sigma;
  }

  inline bool operator != (const CBetaOpts& opts) const {
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

struct DLLEXPORT_QMEAN CBetaSpatialOrganizerItem{

  CBetaSpatialOrganizerItem(ost::conop::AminoAcid a, geom::Vec3 p, int n, const String& c_n):
                       aa(a), pos(p), number(n), chain_name(c_n) { }
  ost::conop::AminoAcid aa;
  geom::Vec3 pos;
  int number;
  String chain_name;

};

class DLLEXPORT_QMEAN CBetaPotentialImpl : public ost::mol::EntityVisitor {
public:

  CBetaPotentialImpl(): env_(5) { }

  CBetaPotentialImpl(CBetaOpts& opts):
    opts_(opts), env_(5) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(ost::conop::AminoAcid a,
                             ost::conop::AminoAcid b,
                             Real dist)=0;

  void SetEnvironment(ost::mol::EntityView& env);

protected:
  
  CBetaOpts opts_;
  ost::mol::SpatialOrganizer<CBetaSpatialOrganizerItem, geom::Vec3>  env_;

};

}} //namespace




#endif // CBETA_IMPL_HH

