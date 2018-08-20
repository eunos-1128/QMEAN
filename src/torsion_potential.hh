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

