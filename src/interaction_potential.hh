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

#ifndef INTERACTION_POTENTIAL_HH
#define INTERACTION_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/interaction_impl.hh>
#include <qmean/interaction_statistics.hh>

namespace qmean {

class InteractionPotential;
typedef boost::shared_ptr<InteractionPotential> InteractionPotentialPtr;
typedef MultiClassifier<float, int, int, float> InteractionEnergies;


class DLLEXPORT_QMEAN InteractionPotential : public PotentialBase, public impl::InteractionPotentialImpl{
public:

  InteractionPotential() { }
  
  static InteractionPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(atom::ChemType a, atom::ChemType b, Real dist);

  static InteractionPotentialPtr Create(InteractionStatisticPtr s, Real sigma, const String& reference_state);

  static InteractionPotentialPtr Create(InteractionStatisticPtr s1, InteractionStatisticPtr s2, Real sigma, const String& reference_state, Real max_energy);

  void Fill(InteractionStatisticPtr stat, const String& reference_state);

  void Fill(InteractionStatisticPtr s1, InteractionStatisticPtr s2, const String& reference_state, Real max_energy);

  PotentialType GetType() const { return Interaction; }

  void SetEnvironment(ost::mol::EntityView& env) { env_=env; }

  Real GetEnergy(atom::ChemType a, atom::ChemType b, Real dist);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() const { return count_; }

  impl::InteractionOpts GetOpts() const { return opts_; }

  InteractionEnergies* Data() { return &energies_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & energies_;
    ds & opts_;
  }

private:
  
  InteractionEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // INTERACTION_POTENTIAL_HH

