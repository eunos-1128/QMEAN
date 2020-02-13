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

#ifndef QMEAN_HBOND_STATISTICS_HH
#define QMEAN_HBOND_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/hbond_impl.hh>

namespace qmean {

class HBondStatistic;
typedef boost::shared_ptr<HBondStatistic> HBondStatisticPtr;
typedef FloatHistogram<int,float,float,float,float> HBondHistogram;
// parametrized as state, dist, alpha, beta, gamma


class DLLEXPORT_QMEAN HBondStatistic : public StatisticBase, impl::HBondPotentialImpl{
public:
  HBondStatistic() { }

  HBondStatistic(Real d_min, Real d_max, Real alpha_min, Real alpha_max,
                 Real beta_min, Real beta_max, Real gamma_min, Real gamma_max, 
                 int d_bins, int alpha_bins, int beta_bins, int gamma_bins,
                 int seq_sep);

  static HBondStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(int state,Real d,Real alpha,Real beta,Real gamma);

  void Extract(ost::mol::EntityView& view, Real weight);

  Real GetTotalCount() const;

  Real GetCount(uint state) const;

  Real GetCount(uint state, uint d_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const;

  PotentialType GetType() const { return HBond; }

  impl::HBondOpts GetOpts() const{ return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  HBondHistogram histo_;
  Real weight_;
};

}//namespace
#endif // HBOND_STATISTICS_HH
