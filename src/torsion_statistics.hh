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

#ifndef TORSION_STATISTICS_HH
#define TORSION_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/torsion_impl.hh>

namespace qmean {

class TorsionStatistic;
typedef boost::shared_ptr<TorsionStatistic> TorsionStatisticPtr;
typedef FloatHistogram<float,float,float,float,float,float> TorsionHistogram;

class DLLEXPORT_QMEAN TorsionStatistic : public StatisticBase, public impl::TorsionPotentialImpl {
public:
  TorsionStatistic() { }

  TorsionStatistic(std::vector<String>& group_identifier, std::vector<int>& num_of_bins);
 
  static TorsionStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(const String& gi, std::vector<Real>& angles);

  void Extract(ost::mol::EntityView& target, Real weight);

  Real GetTotalCount() const;

  Real GetTotalCount(const String& gi) const;

  Real GetCount(const String& gi, std::vector<int>& bins) const;

  Real GetCount(std::vector<int>& bins) const;

  PotentialType GetType() const { return Torsion; }

  impl::TorsionOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;

    if(ds.IsSource()){
      for(std::vector<String>::iterator it=opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){
        histos_[*it]=TorsionHistogram();
      }
    }

    for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
      ds & it->second;
    }
  }

private:

  std::map<String, TorsionHistogram> histos_;
  Real weight_;
};

} //namespace


#endif // TORSION_STATISTICS_HH

