#ifndef QMEAN_PACKING_STATISTICS_HH
#define QMEAN_PACKING_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/packing_impl.hh>

namespace qmean {

class PackingStatistic;
typedef boost::shared_ptr<PackingStatistic> PackingStatisticPtr;
typedef FloatHistogram<int,int> PackingHistogram;

// parametrized as atom, number of atoms in environment


class DLLEXPORT_QMEAN PackingStatistic : public StatisticBase, impl::PackingPotentialImpl{
public:
  PackingStatistic() { }

  PackingStatistic(Real cutoff_radius, int max_count);

  static PackingStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& bins);

  void Extract(ost::mol::EntityView& target, ost::mol::EntityView& env, Real weight);

  Real GetTotalCount() const;

  Real GetCount(atom::ChemType a) const;

  Real GetCount(atom::ChemType a, uint bin) const;

  Real GetCount(uint bin) const;

  PotentialType GetType() const { return Packing; }

  impl::PackingOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & histo_;
    ds & opts_;
  }

private:
  PackingHistogram histo_;
  Real weight_;
};


}//namespace
#endif // PACKING_STATISTICS_HH

