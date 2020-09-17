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

#include <qmean/torsion_statistics.hh>

namespace qmean{

TorsionStatistic::TorsionStatistic(std::vector<String>& gi, std::vector<int>& nob){

  for(size_t i=0;i<gi.size();++i){
    histos_[gi[i]]=TorsionHistogram(ContinuousClassifier(nob[0], -M_PI, M_PI),
                                    ContinuousClassifier(nob[1], -M_PI, M_PI),
                                    ContinuousClassifier(nob[2], -M_PI, M_PI),
                                    ContinuousClassifier(nob[3], -M_PI, M_PI),
                                    ContinuousClassifier(nob[4], -M_PI, M_PI),
                                    ContinuousClassifier(nob[5], -M_PI, M_PI));
  }
  impl::TorsionOpts opts(gi, nob);
  opts_=opts;
}

TorsionStatisticPtr TorsionStatistic::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  TorsionStatisticPtr torsion_p(new TorsionStatistic);
  ds >> *torsion_p;
  return torsion_p;
}

void TorsionStatistic::OnSave(ost::io::BinaryDataSink& ds){
  ds << *this;
}

void TorsionStatistic::OnInteraction(const String& gi, std::vector<Real>& angles){
    histos_[gi].Add(weight_, angles[0], angles[1], angles[2], angles[3], angles[4], angles[5]);
}

void TorsionStatistic::Extract(ost::mol::EntityView& target, Real weight){
  weight_=weight;
  target.Apply(*this);
}

Real TorsionStatistic::GetTotalCount() const {

  Real total=0.0;
  for(std::map<String,TorsionHistogram>::const_iterator it=histos_.begin();it!=histos_.end();++it){
    total+=this->GetTotalCount(it->first);
  }
  return total;
}

Real TorsionStatistic::GetTotalCount(const String& gi) const{

  if(histos_.count(gi)==0){
    std::stringstream ss;
    ss<<"invalid group identifier: ";
    ss<<gi;
    throw std::runtime_error(ss.str());
  }

  typedef TorsionHistogram::IndexType Index;
  Real total=0.0;

  for(int a=0;a<opts_.num_of_bins[0];++a){
    for(int b=0;b<opts_.num_of_bins[1];++b){
      for(int c=0;c<opts_.num_of_bins[2];++c){
        for(int d=0;d<opts_.num_of_bins[3];++d){
          for(int e=0;e<opts_.num_of_bins[4];++e){
            for(int f=0;f<opts_.num_of_bins[5];++f){
              std::map<String,TorsionHistogram>::const_iterator i = histos_.find(gi);
              total += i->second.Get(Index(a,b,c,d,e,f));
            }
          }
        }
      }
    }
  }

  return total;
}

Real TorsionStatistic::GetCount(const String& gi, std::vector<int>& bins) const{

  if(histos_.count(gi)==0){
    std::stringstream ss;
    ss<<"invalid group identifier: ";
    ss<<gi;
    throw std::runtime_error(ss.str());
  }
  if(bins.size()!=6){
    throw std::runtime_error("There are 6 angle bins!");
  }

  for(int i = 0; i < 6; ++i){
    if(bins[i] < 0 || bins[i] >= opts_.num_of_bins[i]){
      throw std::runtime_error("Cannot get count for invalid torsion bin!");
    }
  }
  typedef TorsionHistogram::IndexType Index;
  std::map<String,TorsionHistogram>::const_iterator i = histos_.find(gi);
  return i->second.Get(Index(bins[0],bins[1],bins[2],bins[3],bins[4],bins[5]));
}

Real TorsionStatistic::GetCount(std::vector<int>& bins) const{

  if(bins.size()!=6){
    throw std::runtime_error("There are 6 angle bins!");
  }

  for(int i = 0; i < 6; ++i){
    if(bins[i] < 0 || bins[i] >= opts_.num_of_bins[i]){
      throw std::runtime_error("Cannot get count for invalid torsion bin!");
    }
  }

  Real total=0;

  for(std::map<String,TorsionHistogram>::const_iterator it=histos_.begin();it!=histos_.end();++it){
    total+=this->GetCount(it->first, bins);
  }
  return total;
}

}

