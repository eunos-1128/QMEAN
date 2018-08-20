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

