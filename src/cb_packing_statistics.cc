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

#include <qmean/cb_packing_statistics.hh>

namespace qmean {

CBPackingStatistic::CBPackingStatistic(Real cutoff_radius,
                                       int max_count){
   

  impl::CBPackingOpts opts(cutoff_radius, max_count, 1);
  opts_=opts;
  CBPackingHistogram histo=CBPackingHistogram(IntegralClassifier(ost::conop::XXX, 0),       
                                              IntegralClassifier(max_count+1, 0));
  histo_=histo;
}

CBPackingStatisticPtr CBPackingStatistic::Load(const String& filename){
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  CBPackingStatisticPtr p(new CBPackingStatistic);
  ds >> *p;
  return p;
}

void CBPackingStatistic::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}

void CBPackingStatistic::OnInteraction(ost::conop::AminoAcid aa,
                                       int bin){
  histo_.Add(weight_, aa, bin);
}

void CBPackingStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight)
{
  this->SetEnvironment(env);
  weight_=weight;
  target.Apply(*this);
}

Real CBPackingStatistic::GetTotalCount() const{
  Real total_count=0;
  for(int aa=0;aa<ost::conop::XXX;++aa){
    total_count+=this->GetCount(ost::conop::AminoAcid(aa));
  }
  return total_count;
}

Real CBPackingStatistic::GetCount(ost::conop::AminoAcid aa) const{

  Real total_count=0;
  for(int i=0;i<opts_.max_counts+1;++i){
    total_count += this->GetCount(aa,i);
  }
  return total_count;
}

Real CBPackingStatistic::GetCount(ost::conop::AminoAcid aa, uint bin) const{

  if(aa == ost::conop::XXX){
    throw std::runtime_error("Cannot get count for invalid amino acid!");
  }

  if(bin >= opts_.max_counts+1){
    throw std::runtime_error("Cannot get count for invalid bin!");
  }

  return histo_.Get(CBPackingHistogram::IndexType(aa,bin));
}

Real CBPackingStatistic::GetCount(uint bin) const{

  Real total_count=0;
  for(int i=0; i < ost::conop::XXX; ++i){
    total_count += this->GetCount(ost::conop::AminoAcid(i),bin);
  }
  return total_count;
}

}//namespace
