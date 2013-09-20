#include <qmean/torsion_statistics.hh>

namespace qmean{

TorsionStatistic::TorsionStatistic(std::vector<String>& gi, std::vector<int>& nob){

  for(int i=0;i<gi.size();++i){
    histos_[gi[i]]=TorsionHistogram(ContinuousClassifier(nob[0], -M_PI, M_PI),
                                    ContinuousClassifier(nob[1], -M_PI, M_PI),
                                    ContinuousClassifier(nob[2], -M_PI, M_PI),
                                    ContinuousClassifier(nob[3], -M_PI, M_PI),
                                    ContinuousClassifier(nob[4], -M_PI, M_PI),
                                    ContinuousClassifier(nob[5], -M_PI, M_PI));
  }
  impl::TorsionOpts opts(gi, nob);
  opts_=opts;
  std::cerr<<"finished initialisation..."<<std::endl;
}

TorsionStatisticPtr TorsionStatistic::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  TorsionStatisticPtr torsion_p(new TorsionStatistic);
  ds >> *torsion_p;
  return torsion_p;
}

void TorsionStatistic::OnSave(ost::io::BinaryDataSink& ds){
  ds << *this;
}

void TorsionStatistic::OnInteraction(const String& gi, std::vector<Real>& angles){
    histos_[gi].Add(weight_, angles[0], angles[1], angles[2], angles[3], angles[4], angles[5]);
}

void TorsionStatistic::Extract(ost::mol::EntityView& target, Real weight){
  weight_=weight;
  target.Apply(*this);
}

Real TorsionStatistic::GetTotalCount(){

  Real total=0.0;
  for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
    total+=this->GetTotalCount(it->first);
  }
  return total;
}

Real TorsionStatistic::GetTotalCount(const String& gi){

  if(histos_.count(gi)==0){
    std::stringstream ss;
    ss<<"invalid group identifier: ";
    ss<<gi;
    throw std::runtime_error(ss.str());
  }

  typedef TorsionHistogram::IndexType Index;
  Real total=0.0;

  for(int a=0;a<opts_.num_of_bins[0];++a){
    for(int b=0;b<opts_.num_of_bins[1];++b){
      for(int c=0;c<opts_.num_of_bins[2];++c){
        for(int d=0;d<opts_.num_of_bins[3];++d){
          for(int e=0;e<opts_.num_of_bins[4];++e){
            for(int f=0;f<opts_.num_of_bins[5];++f){
              total+=histos_[gi].Get(Index(a,b,c,d,e,f));
            }
          }
        }
      }
    }
  }

  return total;
}

Real TorsionStatistic::GetCount(const String& gi, std::vector<int>& bins){

  if(histos_.count(gi)==0){
    std::stringstream ss;
    ss<<"invalid group identifier: ";
    ss<<gi;
    for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
      std::cout<<it->first<<std::endl;
    }
    throw std::runtime_error(ss.str());
  }
  if(bins.size()!=6){
    throw std::runtime_error("There are 6 angle bins!");
  }
  typedef TorsionHistogram::IndexType Index;

  return histos_[gi].Get(Index(bins[0],bins[1],bins[2],bins[3],bins[4],bins[5]));
}

Real TorsionStatistic::GetCount(std::vector<int>& bins){

  if(bins.size()!=6){
    throw std::runtime_error("There are 6 angle bins!");
  }

  Real total=0;

  for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
    total+=this->GetCount(it->first, bins);
  }
  return total;
}

}
