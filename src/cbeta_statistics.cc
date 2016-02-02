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

