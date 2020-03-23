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

#ifndef QMEAN_PACKING_STATISTICS_HH
#define QMEAN_PACKING_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/packing_impl.hh>

namespace qmean {

class PackingStatistic;
typedef boost::shared_ptr<PackingStatistic> PackingStatisticPtr;
typedef FloatHistogram<int,int> PackingHistogram;

// parametrized as atom, number of atoms in environment


class DLLEXPORT_QMEAN PackingStatistic : public StatisticBase, impl::PackingPotentialImpl{
public:
  PackingStatistic() { }

  PackingStatistic(Real cutoff_radius, int max_count);

  static PackingStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& bins);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(atom::ChemType a) const;

  Real GetCount(atom::ChemType a, uint bin) const;

  Real GetCount(uint bin) const;

  PotentialType GetType() const { return Packing; }

  impl::PackingOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  PackingHistogram histo_;
  Real weight_;
};


}//namespace
#endif // PACKING_STATISTICS_HH

