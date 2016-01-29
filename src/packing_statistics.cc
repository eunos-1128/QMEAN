#include <qmean/packing_statistics.hh>

namespace qmean {

PackingStatistic::PackingStatistic(Real cutoff_radius,
                                     int max_count){
   
  impl::PackingOpts opts(cutoff_radius, max_count, 1);
  opts_=opts;
  PackingHistogram histo=PackingHistogram(IntegralClassifier(atom::UNKNOWN, 0),       
                                          IntegralClassifier(max_count+1, 0));
  histo_=histo;
}

PackingStatisticPtr PackingStatistic::Load(const String& filename){
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  PackingStatisticPtr p(new PackingStatistic);
  ds >> *p;
  return p;
}

void PackingStatistic::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}

void PackingStatistic::OnInteraction(ost::conop::AminoAcid aa,
                                     const ost::mol::AtomHandleList& a_list,
                                     std::vector<int>& bins){

  ost::mol::AtomHandleList::const_iterator a_it = a_list.begin();
  std::vector<int>::iterator b_it = bins.begin();
  
  for(;a_it!=a_list.end() && b_it!=bins.end();++a_it,++b_it){
    atom::ChemType a = qmean::GetAtomTypeByName(aa, a_it->GetName());
    if(a!=atom::UNKNOWN){
      histo_.Add(weight_, a, *b_it);
    }
  }
}

void PackingStatistic::Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight)
{
  env_=env;
  weight_=weight;
  target.Apply(*this);
}

Real PackingStatistic::GetTotalCount(){
  Real total_count=0;
  for(int a=0;a<atom::UNKNOWN;++a){
    total_count+=this->GetCount(atom::ChemType(a));
  }
  return total_count;
}

Real PackingStatistic::GetCount(atom::ChemType a){

  Real total_count=0;
  for(int i=0;i<opts_.max_counts+1;++i){
    total_count+=histo_.Get(a,i);
  }
  return total_count;
}

Real PackingStatistic::GetCount(atom::ChemType a, int bin){
  return histo_.Get(a,bin);
}

Real PackingStatistic::GetCount(int bin){

  Real total_count=0;
  for(int i=0;i<atom::UNKNOWN;++i){
    total_count+=histo_.Get(i, bin);
  }
  return total_count;
}

}//namespace

