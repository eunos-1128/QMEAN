#ifndef MREDUCED_STATISTICS_HH
#define MREDUCED_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/reduced_impl.hh>

namespace qmean {

class ReducedStatistic;
typedef boost::shared_ptr<ReducedStatistic> ReducedStatisticPtr;
typedef FloatHistogram<int,int,float,float,float,float> ReducedHistogram;

// parametrized as first amino acid, second amino acid,
// calpha-calpha distance, angle


class DLLEXPORT_QMEAN ReducedStatistic : public StatisticBase, impl::ReducedPotentialImpl{
public:
  ReducedStatistic() { }

  ReducedStatistic(Real l_cutoff, Real u_cutoff, int nab, int ndab,
                    int ndb, int ssep);

  static ReducedStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta, Real gamma);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(ost::conop::AminoAcid aa_one,
                ost::conop::AminoAcid aa_two) const;


  Real GetCount(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const;

  Real GetCount(uint dist_bin, uint alpha_bin, uint beta_bin, uint gamma_bin) const;

  PotentialType GetType() const { return Reduced; }

  impl::ReducedOpts GetOpts() const { return opts_; }

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

