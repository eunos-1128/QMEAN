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

#include <qmean/reduced_statistics.hh>

namespace qmean {

ReducedStatistic::ReducedStatistic(Real l_cutoff, Real u_cutoff, int nab,
                                   int ndab, int ndb, int ssep):
    histo_(IntegralClassifier(20,0), 
           IntegralClassifier(20,0), 
           ContinuousClassifier(ndb, l_cutoff, u_cutoff), 
           ContinuousClassifier(nab, 0.0, M_PI),
           ContinuousClassifier(nab, 0.0, M_PI),
           ContinuousClassifier(ndab, -M_PI, M_PI)){

    impl::ReducedOpts opts(l_cutoff, u_cutoff, nab, ndab, ndb, ssep);
    opts_=opts;
}

ReducedStatisticPtr ReducedStatistic::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  ReducedStatisticPtr reduced_p(new ReducedStatistic);
  ds >> *reduced_p;
  return reduced_p;
}

void ReducedStatistic::OnSave(ost::io::BinaryDataSink& ds){
  ds << *this;
}

void ReducedStatistic::OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                                     Real dist, Real alpha, Real beta, Real gamma){ 
  histo_.Add(weight_, aa_one, aa_two, dist, alpha, beta, gamma);
}

void ReducedStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight){
  weight_=weight;
  this->SetEnvironment(env);
  target.Apply(*this);
}

Real ReducedStatistic::GetTotalCount() const{

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (int i=0; i<ost::conop::XXX; ++i) {
    for (int j=0; j<ost::conop::XXX; ++j) {
      for (int k=0; k<opts_.num_dist_bins; ++k) {
        for (int l=0; l<opts_.num_angle_bins; ++l) {
          for (int m=0; m<opts_.num_angle_bins; ++m) {
            for (int n=0; n<opts_.num_dihedral_bins; ++n) {
              count+=histo_.Get(Index(i, j, k, l, m, n));
            }
          }
        }
      }
    }
  }
  return count;
}

Real ReducedStatistic::GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two) const{

  if(aa_one == ost::conop::XXX || aa_two == ost::conop::XXX){
    throw std::runtime_error("Cannot get count of invalid amino acid!"); 
  }

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (int k=0; k<opts_.num_dist_bins; ++k) {
    for (int l=0; l<opts_.num_angle_bins; ++l) {
      for (int m=0; m<opts_.num_angle_bins; ++m) {
        for (int n=0; n<opts_.num_dihedral_bins; ++n) {
          count+=histo_.Get(Index(aa_one, aa_two, k, l, m, n));
        }
      }
    }
  }
  return count;
}

Real ReducedStatistic::GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                                uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const{

  if(aa_one == ost::conop::XXX || aa_two == ost::conop::XXX){
    throw std::runtime_error("Cannot get count of invalid amino acid!"); 
  }

  if(static_cast<int>(dist_bin) >= opts_.num_dist_bins){
    throw std::runtime_error("Cannot get count of invalid distance bin!");
  }

  if(static_cast<int>(alpha_bin) >= opts_.num_angle_bins){
    throw std::runtime_error("Cannot get count of invalid alpha bin!");
  }

  if(static_cast<int>(beta_bin) >= opts_.num_angle_bins){
    throw std::runtime_error("Cannot get count of invalid beta bin!");
  }

  if(static_cast<int>(gamma_bin) >= opts_.num_dihedral_bins){
    throw std::runtime_error("Cannot get count of invalid gamma bin!");
  }

  return histo_.Get(ReducedHistogram::IndexType(aa_one, aa_two, dist_bin,
                                                alpha_bin, beta_bin, gamma_bin));
}

Real ReducedStatistic::GetCount(uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const{

  if(static_cast<int>(dist_bin) >= opts_.num_dist_bins){
    throw std::runtime_error("Cannot get count of invalid distance bin!");
  }

  if(static_cast<int>(alpha_bin) >= opts_.num_angle_bins){
    throw std::runtime_error("Cannot get count of invalid alpha bin!");
  }

  if(static_cast<int>(beta_bin) >= opts_.num_angle_bins){
    throw std::runtime_error("Cannot get count of invalid beta bin!");
  }

  if(static_cast<int>(gamma_bin) >= opts_.num_dihedral_bins){
    throw std::runtime_error("Cannot get count of invalid gamma bin!");
  }

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (size_t i=0; i<ost::conop::XXX; ++i) {
    for (size_t j=0; j<ost::conop::XXX; ++j) {
      count+=histo_.Get(Index(i, j, dist_bin, alpha_bin, beta_bin, gamma_bin));
    }
  }
  return count;
}

}

