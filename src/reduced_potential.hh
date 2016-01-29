#ifndef MREDUCED_POTENTIAL_HH
#define MREDUCED_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/reduced_impl.hh>
#include <qmean/reduced_statistics.hh>

namespace qmean {

class ReducedPotential;
typedef boost::shared_ptr<ReducedPotential> ReducedPotentialPtr;
typedef MultiClassifier<float, int, int, float, float, float, float> ReducedEnergies;


class DLLEXPORT_QMEAN ReducedPotential : public PotentialBase, public impl::ReducedPotentialImpl{
public:

  ReducedPotential() { }
  
  static ReducedPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta, Real gamma);

  static ReducedPotentialPtr Create(ReducedStatisticPtr s, Real sigma, const String& reference_state);

  void Fill(ReducedStatisticPtr stat, const String& reference_state);

  PotentialType GetType() { return Reduced; }

  Real GetEnergy(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real distance, Real alpha, Real beta, Real gamma);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() { return count_; }

  impl::ReducedOpts GetOpts() { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  ReducedEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // MREDUCED_POTENTIAL_HH

