#include <qmean/reduced_statistics.hh>

namespace qmean {

ReducedStatistic::ReducedStatistic(Real l_cutoff, Real u_cutoff, int nab,
                                   int ndb, int ssep):
    histo_(IntegralClassifier(20,0), IntegralClassifier(20,0), ContinuousClassifier(ndb, l_cutoff, u_cutoff), ContinuousClassifier(nab, 0.0, M_PI)){

    impl::ReducedOpts opts(l_cutoff, u_cutoff, nab, ndb, ssep);
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
                                                 Real dist, Real angle){ 
  histo_.Add(weight_, aa_one, aa_two, dist, angle);
}

void ReducedStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight){
  weight_=weight;
  this->SetEnvironment(env);
  target.Apply(*this);
}

Real ReducedStatistic::GetTotalCount(){

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (size_t i=0; i<ost::conop::XXX; ++i) {
    for (size_t j=0; j<ost::conop::XXX; ++j) {
      for (size_t k=0; k<opts_.num_dist_bins; ++k) {
        for (size_t l=0; l<opts_.num_angular_bins; ++l) {
          count+=histo_.Get(Index(i, j, k, l));
        }
      }
    }
  }
  return count;
}

Real ReducedStatistic::GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two){

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (size_t k=0; k<opts_.num_dist_bins; ++k) {
    for (size_t l=0; l<opts_.num_angular_bins; ++l) {
      count+=histo_.Get(Index(aa_one, aa_two, k, l));
    }
  }
  return count;
}

Real ReducedStatistic::GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
              int dist_bin, int ang_bin){

  return histo_.Get(ReducedHistogram::IndexType(aa_one, aa_two, dist_bin,
                                                ang_bin));
}

Real ReducedStatistic::GetCount(int dist_bin, int ang_bin) {

  typedef ReducedHistogram::IndexType Index;
  Real count=0.0;
  for (size_t i=0; i<ost::conop::XXX; ++i) {
    for (size_t j=0; j<ost::conop::XXX; ++j) {
      count+=histo_.Get(Index(i, j, dist_bin, ang_bin));
    }
  }
  return count;
}

}

