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

#ifndef QMEAN_PACKING_POTENTIAL_HH
#define QMEAN_PACKING_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/packing_impl.hh>
#include <qmean/packing_statistics.hh>

namespace qmean {

class PackingPotential;
typedef boost::shared_ptr<PackingPotential> PackingPotentialPtr;
typedef MultiClassifier<float, int, int> PackingEnergies;


class DLLEXPORT_QMEAN PackingPotential : public PotentialBase, public impl::PackingPotentialImpl{
public:

  PackingPotential() { }
  
  static PackingPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& counts);

  static PackingPotentialPtr Create(PackingStatisticPtr stat, Real sigma, const String& reference_state);

  void Fill(PackingStatisticPtr stat, const String& reference_state);

  PotentialType GetType() const { return Packing; }

  void SetEnvironment(ost::mol::EntityView& env) { env_=env; }

  Real GetEnergy(atom::ChemType a, int count);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  impl::PackingOpts GetOpts() const { return opts_; }

  PackingEnergies* Data() { return &energies_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  PackingEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // PACKING_POTENTIAL_HH

