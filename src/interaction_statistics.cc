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

#include <qmean/interaction_statistics.hh>

namespace qmean {

InteractionStatistic::InteractionStatistic(Real l_cutoff,
                                             Real u_cutoff,
                                             int nob,
                                             int ssep):

                   histo_(IntegralClassifier(atom::UNKNOWN, 0),
                          IntegralClassifier(atom::UNKNOWN, 0),
                          ContinuousClassifier(nob, l_cutoff, u_cutoff))
{
  impl::InteractionOpts opts(l_cutoff, u_cutoff, nob, ssep);
  opts_=opts;
}

InteractionStatisticPtr InteractionStatistic::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  InteractionStatisticPtr p(new InteractionStatistic);
  ds >> *p;
  return p;
}

void InteractionStatistic::OnSave(ost::io::BinaryDataSink& ds){
  ds << *this;
}

void InteractionStatistic::OnInteraction(atom::ChemType a, atom::ChemType b, Real dist){
  histo_.Add(weight_, a, b, dist);
}

void InteractionStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight){
  env_=env;
  weight_=weight;
  target.Apply(*this);
}

Real InteractionStatistic::GetTotalCount() const{
  Real total_count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    total_count+=this->GetCount(i);
  }
  return total_count;
}

Real InteractionStatistic::GetCount(atom::ChemType a, atom::ChemType b, uint dist_bin) const{

  if(a == atom::UNKNOWN || b == atom::UNKNOWN){
    throw std::runtime_error("Cannot get count for invalid ChemType!");
  }

  if(dist_bin >= opts_.number_of_bins){
    throw std::runtime_error("Cannot get count for invalid distance bin!");
  }

  return histo_.Get(InteractionHistogram::IndexType(a, b, dist_bin)); 
}

Real InteractionStatistic::GetCount(atom::ChemType a, atom::ChemType b) const{
  Real count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    count+=this->GetCount(a,b,i);
  }
  return count;
}

Real InteractionStatistic::GetCount(uint dist_bin) const{

  Real count=0.0;
  for (size_t i=0; i<atom::UNKNOWN; ++i) {
    for (size_t j=0; j<atom::UNKNOWN; ++j) {
      count+=this->GetCount(atom::ChemType(i),atom::ChemType(j),dist_bin);
    }
  }
  return count;
}

}

