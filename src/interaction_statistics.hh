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

#ifndef INTERACTION_STATISTICS_HH
#define INTERACTION_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/interaction_impl.hh>

namespace qmean {

class InteractionStatistic;
typedef boost::shared_ptr<InteractionStatistic> InteractionStatisticPtr;
typedef FloatHistogram<int, int, float> InteractionHistogram;

// parametrized as first atom, second atom, distance

class DLLEXPORT_QMEAN InteractionStatistic : public StatisticBase, impl::InteractionPotentialImpl{
public:
  InteractionStatistic() { }

  InteractionStatistic(Real l_cutoff, Real u_cutoff, int nob, int ssep);

  static InteractionStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(atom::ChemType a, atom::ChemType b, Real dist);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(atom::ChemType a, atom::ChemType b, uint dist_bin) const;

  Real GetCount(atom::ChemType a, atom::ChemType b) const;

  Real GetCount(uint dist_bin) const;

  PotentialType GetType() const { return Interaction; }

  impl::InteractionOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  InteractionHistogram histo_;
  Real weight_;
};


}//namespace
#endif // INTERACTION_STATISTICS_HH

