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

#include <qmean/cbeta_statistics.hh>

namespace qmean {

CBetaStatistic::CBetaStatistic(Real l_cutoff,
                               Real u_cutoff,
                               int nob,
                               int ssep):

                   histo_(IntegralClassifier(ost::conop::XXX, 0),
                          IntegralClassifier(ost::conop::XXX, 0),
                          ContinuousClassifier(nob, l_cutoff, u_cutoff))
{
  impl::CBetaOpts opts(l_cutoff, u_cutoff, nob, ssep);
  opts_=opts;
}

CBetaStatisticPtr CBetaStatistic::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  CBetaStatisticPtr p(new CBetaStatistic);
  ds >> *p;
  return p;
}

void CBetaStatistic::OnSave(ost::io::BinaryDataSink& ds){
  ds << *this;
}

void CBetaStatistic::OnInteraction(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist){
  histo_.Add(weight_, a, b, dist);
}

void CBetaStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight){
  this->SetEnvironment(env);
  weight_=weight;
  target.Apply(*this);
}

Real CBetaStatistic::GetTotalCount() const{
  Real total_count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    total_count+=this->GetCount(i);
  }
  return total_count;
}

Real CBetaStatistic::GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b, uint dist_bin) const{

  if(a == ost::conop::XXX || b == ost::conop::XXX){
    throw std::runtime_error("Cannot get count for invalid AminoAcid!");
  }

  if(dist_bin >= opts_.number_of_bins){
    throw std::runtime_error("Cannot get count for invalid bin!");
  }

  return histo_.Get(CBetaHistogram::IndexType(a, b, dist_bin)); 
}

Real CBetaStatistic::GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b) const{

  Real count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    count+=this->GetCount(a,b,i);
  }
  return count;
}

Real CBetaStatistic::GetCount(uint dist_bin) const{

  Real count=0.0;
  for (size_t i=0; i<ost::conop::XXX; ++i) {
    for (size_t j=0; j<ost::conop::XXX; ++j) {
      count+=this->GetCount(ost::conop::AminoAcid(i),ost::conop::AminoAcid(j),dist_bin);
    }
  }
  return count;
}

}

