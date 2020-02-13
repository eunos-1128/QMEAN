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

#ifndef CBETA_POTENTIAL_HH
#define CBETA_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/cbeta_impl.hh>
#include <qmean/cbeta_statistics.hh>


namespace qmean {

class CBetaPotential;
typedef boost::shared_ptr<CBetaPotential> CBetaPotentialPtr;
typedef MultiClassifier<float, int, int, float> CBetaEnergies;


class DLLEXPORT_QMEAN CBetaPotential : public PotentialBase, public impl::CBetaPotentialImpl{
public:

  CBetaPotential() { }
  
  static CBetaPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist);

  static CBetaPotentialPtr Create(CBetaStatisticPtr s, Real sigma, const String& reference_state);

  static CBetaPotentialPtr Create(CBetaStatisticPtr s1, CBetaStatisticPtr s2, Real sigma, const String& reference_state, Real max_energy);

  void Fill(CBetaStatisticPtr stat, const String& reference_state);

  void Fill(CBetaStatisticPtr s1, CBetaStatisticPtr s2, const String& reference_state, Real max_energy);

  PotentialType GetType() const { return CBeta; }

  Real GetEnergy(ost::conop::AminoAcid a, ost::conop::AminoAcid, Real dist);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() const { return count_; }

  impl::CBetaOpts GetOpts() const { return opts_; }

  CBetaEnergies* Data() { return &energies_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & energies_;
    ds & opts_;
  }

private:


  
  CBetaEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // CBETA_POTENTIAL_HH

