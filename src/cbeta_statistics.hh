// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
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

#ifndef CBETA_STATISTICS_HH
#define CBETA_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/cbeta_impl.hh>

namespace qmean {

class CBetaStatistic;
typedef boost::shared_ptr<CBetaStatistic> CBetaStatisticPtr;
typedef FloatHistogram<int, int, float> CBetaHistogram;

// parametrized as first atom, second atom, distance

class DLLEXPORT_QMEAN CBetaStatistic : public StatisticBase, impl::CBetaPotentialImpl{
public:
  CBetaStatistic() { }

  CBetaStatistic(Real l_cutoff, Real u_cutoff, int nob, int ssep);

  static CBetaStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b, uint dist_bin) const;

  Real GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b) const;

  Real GetCount(uint dist_bin) const;

  PotentialType GetType() const { return CBeta; }

  impl::CBetaOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & histo_;
  }

private:
  CBetaHistogram histo_;
  Real weight_;
};


}//namespace
#endif // CBETA_STATISTICS_HH

