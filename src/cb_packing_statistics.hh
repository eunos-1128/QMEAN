#ifndef QMEAN_CB_PACKING_STATISTICS_HH
#define QMEAN_CB_PACKING_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/cb_packing_impl.hh>

namespace qmean {

class CBPackingStatistic;
typedef boost::shared_ptr<CBPackingStatistic> CBPackingStatisticPtr;
typedef FloatHistogram<int,int> CBPackingHistogram;

// parametrized as residue, number of cb atoms in environment


class DLLEXPORT_QMEAN CBPackingStatistic : public StatisticBase, impl::CBPackingPotentialImpl{
public:
  CBPackingStatistic() { }

  CBPackingStatistic(Real cutoff_radius, int max_count);

  static CBPackingStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             int bin);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(ost::conop::AminoAcid aa) const;

  Real GetCount(ost::conop::AminoAcid aa, uint bin) const;

  Real GetCount(uint bin) const;

  PotentialType GetType() const { return CBPacking; }

  impl::CBPackingOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  CBPackingHistogram histo_;
  Real weight_;
};


}//namespace
#endif // CB_PACKING_STATISTICS_HH

