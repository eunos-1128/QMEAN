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

#ifndef QMEAN_CB_PACKING_STATISTICS_HH
#define QMEAN_CB_PACKING_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/cb_packing_impl.hh>

namespace qmean {

class CBPackingStatistic;
typedef boost::shared_ptr<CBPackingStatistic> CBPackingStatisticPtr;
typedef FloatHistogram<int,int> CBPackingHistogram;

// parametrized as residue, number of cb atoms in environment


class DLLEXPORT_QMEAN CBPackingStatistic : public StatisticBase, impl::CBPackingPotentialImpl{
public:
  CBPackingStatistic() { }

  CBPackingStatistic(Real cutoff_radius, int max_count);

  static CBPackingStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             int bin);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(ost::conop::AminoAcid aa) const;

  Real GetCount(ost::conop::AminoAcid aa, uint bin) const;

  Real GetCount(uint bin) const;

  PotentialType GetType() const { return CBPacking; }

  impl::CBPackingOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  CBPackingHistogram histo_;
  Real weight_;
};


}//namespace
#endif // CB_PACKING_STATISTICS_HH

