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

