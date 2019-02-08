// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
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

