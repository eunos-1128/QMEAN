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

