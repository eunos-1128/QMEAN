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

#ifndef MREDUCED_STATISTICS_HH
#define MREDUCED_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/reduced_impl.hh>

namespace qmean {

class ReducedStatistic;
typedef boost::shared_ptr<ReducedStatistic> ReducedStatisticPtr;
typedef FloatHistogram<int,int,float,float,float,float> ReducedHistogram;

// parametrized as first amino acid, second amino acid,
// calpha-calpha distance, angle


class DLLEXPORT_QMEAN ReducedStatistic : public StatisticBase, impl::ReducedPotentialImpl{
public:
  ReducedStatistic() { }

  ReducedStatistic(Real l_cutoff, Real u_cutoff, int nab, int ndab,
                    int ndb, int ssep);

  static ReducedStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta, Real gamma);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(ost::conop::AminoAcid aa_one,
                ost::conop::AminoAcid aa_two) const;


  Real GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const;

  Real GetCount(uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const;

  PotentialType GetType() const { return Reduced; }

  impl::ReducedOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & histo_;
  }

private:
  ReducedHistogram histo_;
  Real weight_;
};


}//namespace
#endif // MREDUCED_STATISTICS_HH

