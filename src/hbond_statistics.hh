#ifndef QMEAN_HBOND_STATISTICS_HH
#define QMEAN_HBOND_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/hbond_impl.hh>

namespace qmean {

class HBondStatistic;
typedef boost::shared_ptr<HBondStatistic> HBondStatisticPtr;
typedef FloatHistogram<int,float,float,float,float> HBondHistogram;
// parametrized as state, dist, alpha, beta, gamma


class DLLEXPORT_QMEAN HBondStatistic : public StatisticBase, impl::HBondPotentialImpl{
public:
  HBondStatistic() { }

  HBondStatistic(Real d_min, Real d_max, Real alpha_min, Real alpha_max,
                 Real beta_min, Real beta_max, Real gamma_min, Real gamma_max, 
                 int d_bins, int alpha_bins, int beta_bins, int gamma_bins,
                 int seq_sep);

  static HBondStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(int state,Real d,Real alpha,Real beta,Real gamma);

  void Extract(ost::mol::EntityView& view, Real weight);

  Real GetTotalCount();

  Real GetCount(int state);

  Real GetCount(int state, int d_bin, int alpha_bin, int beta_bin, int gamma_bin);

  PotentialType GetType() { return HBond; }

  impl::HBondOpts GetOpts() { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  HBondHistogram histo_;
  Real weight_;
};

}//namespace
#endif // HBOND_STATISTICS_HH
