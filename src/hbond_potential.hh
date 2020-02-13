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

#ifndef QMEAN_HBOND_POTENTIAL_HH
#define QMEAN_HBOND_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/hbond_impl.hh>
#include <qmean/hbond_statistics.hh>

namespace qmean {

class HBondPotential;
typedef boost::shared_ptr<HBondPotential> HBondPotentialPtr;
typedef MultiClassifier<float, int, float, float, float, float> HBondEnergies;


class DLLEXPORT_QMEAN HBondPotential : public PotentialBase, public impl::HBondPotentialImpl{
public:

  HBondPotential() { }
  
  static HBondPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(int state, Real d, Real alpha, 
                             Real beta, Real gamma);

  static HBondPotentialPtr Create(HBondStatisticPtr stat);

  void Fill(HBondStatisticPtr stat);

  PotentialType GetType() const { return HBond; }

  Real GetEnergy(int state, float d, float alpha, float beta, float gamma);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env);

  Real GetEnergy(ost::mol::ResidueView& target);

  std::vector<Real> GetEnergies(ost::mol::EntityView& view);

  Real GetTotalEnergy(ost::mol::EntityView& view, bool normalize);

  impl::HBondOpts GetOpts() const { return opts_; }

  HBondEnergies* Data() { return &energies_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  HBondEnergies energies_;
  Real energy_;
};

}//namespace
#endif // HBOND_POTENTIAL
