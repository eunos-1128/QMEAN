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
