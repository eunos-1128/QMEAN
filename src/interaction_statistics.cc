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

Real InteractionStatistic::GetTotalCount(){
  Real total_count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    total_count+=this->GetCount(i);
  }
  return total_count;
}

Real InteractionStatistic::GetCount(atom::ChemType a, atom::ChemType b, int dist_bin){
  return histo_.Get(InteractionHistogram::IndexType(a, b, dist_bin)); 
}

Real InteractionStatistic::GetCount(atom::ChemType a, atom::ChemType b){
  Real count=0.0;
  for(int i=0; i<opts_.number_of_bins; ++i){
    count+=this->GetCount(a,b,i);
  }
  return count;
}

Real InteractionStatistic::GetCount(int dist_bin){
  Real count=0.0;
  for (size_t i=0; i<atom::UNKNOWN; ++i) {
    for (size_t j=0; j<atom::UNKNOWN; ++j) {
      count+=this->GetCount(atom::ChemType(i),atom::ChemType(j),dist_bin);
    }
  }
  return count;
}

}
