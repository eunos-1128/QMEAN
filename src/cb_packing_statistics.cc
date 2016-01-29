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

Real CBPackingStatistic::GetTotalCount(){
  Real total_count=0;
  for(int aa=0;aa<ost::conop::XXX;++aa){
    total_count+=this->GetCount(ost::conop::AminoAcid(aa));
  }
  return total_count;
}

Real CBPackingStatistic::GetCount(ost::conop::AminoAcid aa){

  Real total_count=0;
  for(int i=0;i<opts_.max_counts+1;++i){
    total_count += this->GetCount(aa,i);
  }
  return total_count;
}

Real CBPackingStatistic::GetCount(ost::conop::AminoAcid aa, int bin){
  return histo_.Get(CBPackingHistogram::IndexType(aa,bin));
}

Real CBPackingStatistic::GetCount(int bin){

  Real total_count=0;
  for(int i=0; i < ost::conop::XXX; ++i){
    total_count += this->GetCount(ost::conop::AminoAcid(i),bin);
  }
  return total_count;
}

}//namespace
