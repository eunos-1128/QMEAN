// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
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

#ifndef TORSION_POTENTIAL_HH
#define TORSION_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/torsion_impl.hh>
#include <qmean/torsion_statistics.hh>

namespace qmean {

class TorsionPotential;
typedef boost::shared_ptr<TorsionPotential> TorsionPotentialPtr;
typedef MultiClassifier<float, float, float, float, float, float, float> TorsionEnergies;

class DLLEXPORT_QMEAN TorsionPotential : public PotentialBase, public impl::TorsionPotentialImpl{

public:

  TorsionPotential() {}

  static TorsionPotentialPtr Create(TorsionStatisticPtr& stat, Real sigma, const String& reference_state);

  static TorsionPotentialPtr Load(const String& filename);

  void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(const String& group_identifier, std::vector<Real>& angles);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target);

  Real GetEnergy(ost::mol::ResidueView& res);

  Real GetEnergy(std::vector<String>& residue_names, std::vector<Real>& angles);

  Real GetEnergy(const String& group_id, std::vector<Real>& angles);

  Real GetTotalEnergy(ost::mol::EntityView& target, bool normalize);

  PotentialType GetType() const { return Torsion; }

  impl::TorsionOpts GetOpts() const { return opts_; }

  TorsionEnergies* Data(std::vector<String>& residue_names);

  template <typename DS>
  void Serialize(DS& ds){

    ds & opts_;

    if(ds.IsSource()){
      for(std::vector<String>::iterator it= opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){
        energies_[*it]=TorsionEnergies();
      }
    }

    for(std::map<String, TorsionEnergies>::iterator it=energies_.begin();it!=energies_.end();++it){
      ds & it->second;
    }
  }

private:

  void Fill(TorsionStatisticPtr stat, const String& reference_state);

  std::map<String, TorsionEnergies> energies_;
  Real energy_;
  uint32_t count_;

};

}


#endif // TORSION_POTENTIAL_HH

