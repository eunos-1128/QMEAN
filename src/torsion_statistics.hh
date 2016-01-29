#ifndef TORSION_STATISTICS_HH
#define TORSION_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/torsion_impl.hh>

namespace qmean {

class TorsionStatistic;
typedef boost::shared_ptr<TorsionStatistic> TorsionStatisticPtr;
typedef FloatHistogram<float,float,float,float,float,float> TorsionHistogram;

class DLLEXPORT_QMEAN TorsionStatistic : public StatisticBase, public impl::TorsionPotentialImpl {
public:
  TorsionStatistic() { }

  TorsionStatistic(std::vector<String>& group_identifier, std::vector<int>& num_of_bins);
 
  static TorsionStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(const String& gi, std::vector<Real>& angles);

  void Extract(ost::mol::EntityView& target, Real weight);

  Real GetTotalCount();

  Real GetTotalCount(const String& gi);

  Real GetCount(const String& gi, std::vector<int>& bins);

  Real GetCount(std::vector<int>& bins);

  PotentialType GetType() { return Torsion; }

  impl::TorsionOpts GetOpts() { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;

    if(ds.IsSource()){
      for(std::vector<String>::iterator it=opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){
        histos_[*it]=TorsionHistogram();
      }
    }

    for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
      ds & it->second;
    }
  }

private:

  std::map<String, TorsionHistogram> histos_;
  Real weight_;
};

} //namespace


#endif // TORSION_STATISTICS_HH

