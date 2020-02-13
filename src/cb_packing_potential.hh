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

#ifndef QMEAN_CB_PACKING_POTENTIAL_HH
#define QMEAN_CB_PACKING_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/cb_packing_impl.hh>
#include <qmean/cb_packing_statistics.hh>

namespace qmean {

class CBPackingPotential;
typedef boost::shared_ptr<CBPackingPotential> CBPackingPotentialPtr;
typedef MultiClassifier<float, int, int> CBPackingEnergies;


class DLLEXPORT_QMEAN CBPackingPotential : public PotentialBase, public impl::CBPackingPotentialImpl{
public:

  CBPackingPotential() { }
  
  static CBPackingPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             int bin);

  static CBPackingPotentialPtr Create(CBPackingStatisticPtr stat, Real sigma, const String& reference_state);

  void Fill(CBPackingStatisticPtr stat, const String& reference_state);

  void Fill(CBPackingStatisticPtr s1, CBPackingStatisticPtr s2, const String& reference_state, Real max_energy);

  PotentialType GetType() const { return CBPacking; }

  Real GetEnergy(ost::conop::AminoAcid aa, int count);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env);

  Real GetEnergy(ost::mol::ResidueView& target);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  impl::CBPackingOpts GetOpts() const { return opts_; }

  CBPackingEnergies* Data() { return &energies_; } 

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  CBPackingEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // CB_PACKING_POTENTIAL_HH

