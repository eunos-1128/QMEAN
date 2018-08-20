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

#ifndef MREDUCED_POTENTIAL_HH
#define MREDUCED_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/reduced_impl.hh>
#include <qmean/reduced_statistics.hh>

namespace qmean {

class ReducedPotential;
typedef boost::shared_ptr<ReducedPotential> ReducedPotentialPtr;
typedef MultiClassifier<float, int, int, float, float, float, float> ReducedEnergies;


class DLLEXPORT_QMEAN ReducedPotential : public PotentialBase, public impl::ReducedPotentialImpl{
public:

  ReducedPotential() { }
  
  static ReducedPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta, Real gamma);

  static ReducedPotentialPtr Create(ReducedStatisticPtr s, Real sigma, const String& reference_state);

  void Fill(ReducedStatisticPtr stat, const String& reference_state);

  PotentialType GetType() const { return Reduced; }

  Real GetEnergy(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real distance, Real alpha, Real beta, Real gamma);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() const { return count_; }

  impl::ReducedOpts GetOpts() const { return opts_; }

  ReducedEnergies* Data() { return &energies_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  ReducedEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // MREDUCED_POTENTIAL_HH

