#ifndef CBETA_STATISTICS_HH
#define CBETA_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/cbeta_impl.hh>

namespace qmean {

class CBetaStatistic;
typedef boost::shared_ptr<CBetaStatistic> CBetaStatisticPtr;
typedef FloatHistogram<int, int, float> CBetaHistogram;

// parametrized as first atom, second atom, distance

class DLLEXPORT_QMEAN CBetaStatistic : public StatisticBase, impl::CBetaPotentialImpl{
public:
  CBetaStatistic() { }

  CBetaStatistic(Real l_cutoff, Real u_cutoff, int nob, int ssep);

  static CBetaStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount();

  Real GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b, int dist_bin);

  Real GetCount(ost::conop::AminoAcid a, ost::conop::AminoAcid b);

  Real GetCount(int dist_bin);

  impl::CBetaOpts& GetOpts() { return opts_; }

  PotentialType GetType() { return CBeta; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & histo_;
  }

private:
  CBetaHistogram histo_;
  Real weight_;
};


}//namespace
#endif // CBETA_STATISTICS_HH