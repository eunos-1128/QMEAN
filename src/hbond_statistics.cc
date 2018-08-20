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

#include <qmean/hbond_statistics.hh>

namespace qmean {

HBondStatistic::HBondStatistic(Real d_min, Real d_max, Real alpha_min, Real alpha_max,
                               Real beta_min, Real beta_max, Real gamma_min, Real gamma_max, 
                               int d_bins, int alpha_bins, int beta_bins, int gamma_bins,
                               int seq_sep): histo_(IntegralClassifier(3, 0),       
                                                    ContinuousClassifier(d_bins, d_min, d_max),
                                                    ContinuousClassifier(alpha_bins, alpha_min, alpha_max),
                                                    ContinuousClassifier(beta_bins, beta_min, beta_max),
                                                    ContinuousClassifier(gamma_bins, gamma_min, gamma_max)) { 

  impl::HBondOpts opts(d_min, d_max, alpha_min, alpha_max, beta_min, beta_max,
                       gamma_min, gamma_max, d_bins, alpha_bins, beta_bins,
                       gamma_bins, seq_sep);
  opts_ = opts;
}

HBondStatisticPtr HBondStatistic::Load(const String& filename){
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  HBondStatisticPtr p(new HBondStatistic);
  ds >> *p;
  return p;
}

void HBondStatistic::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}

void HBondStatistic::OnInteraction(int state, Real d, Real alpha, Real beta, Real gamma){
  histo_.Add(weight_, state, d, alpha, beta, gamma);
}

void HBondStatistic::Extract(ost::mol::EntityView& view, Real weight){
  this->SetEnvironment(view);
  weight_ = weight;
  unique_bonds_ = true;
  view.Apply(*this);
}

Real HBondStatistic::GetTotalCount() const{
  Real total_count=0;
  for(int i=0; i < 3; ++i){
    total_count+=this->GetCount(i);
  }
  return total_count;
}

Real HBondStatistic::GetCount(uint state) const{

  Real total_count=0;
  for(int i = 0; i < opts_.d_bins; ++i){
    for(int j = 0; j < opts_.alpha_bins; ++j){
      for(int k = 0; k < opts_.beta_bins; ++k){
        for(int l = 0; l < opts_.gamma_bins; ++l){
          total_count += this->GetCount(state,i,j,k,l);
        }
      }
    }
  }
  return total_count;
}

Real HBondStatistic::GetCount(uint state, uint d_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const{

  if(state > 2){
    throw std::runtime_error("Cannot get count for invalid state!");
  }

  if(alpha_bin >= opts_.alpha_bins){
    throw std::runtime_error("Cannot get count for invalid alpha bin!");
  }

  if(beta_bin >= opts_.beta_bins){
    throw std::runtime_error("Cannot get count for invalid beta bin!");
  }

  if(gamma_bin >= opts_.gamma_bins){
    throw std::runtime_error("Cannot get count for invalid gamma bin!");
  }

  return histo_.Get(HBondHistogram::IndexType(state,d_bin,alpha_bin,beta_bin,gamma_bin));
}

}//namespace
