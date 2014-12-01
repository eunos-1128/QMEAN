#ifndef MREDUCED_STATISTICS_HH
#define MREDUCED_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/reduced_impl.hh>

namespace qmean {

class ReducedStatistic;
typedef boost::shared_ptr<ReducedStatistic> ReducedStatisticPtr;
typedef FloatHistogram<int,int,float,float> ReducedHistogram;

// parametrized as first amino acid, second amino acid,
// calpha-calpha distance, angle


class DLLEXPORT_QMEAN ReducedStatistic : public StatisticBase, impl::ReducedPotentialImpl{
public:
  ReducedStatistic() { }

  ReducedStatistic(Real l_cutoff, Real u_cutoff, int nab,
                    int ndb, int ssep);

  static ReducedStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real angle);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount();

  Real GetCount(ost::conop::AminoAcid aa_one,
                    ost::conop::AminoAcid aa_two);


  Real GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                int dist_bin, int ang_bin);

  Real GetCount(int dist_bin, int ang_bin);

  impl::ReducedOpts& GetOpts() { return opts_; }

  PotentialType GetType() { return Reduced; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & histo_;
  }

private:
  ReducedHistogram histo_;
  Real weight_;
};


}//namespace
#endif // MREDUCED_STATISTICS_HH

