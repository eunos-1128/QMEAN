#ifndef QMEAN_PACKING_POTENTIAL_HH
#define QMEAN_PACKING_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/packing_impl.hh>
#include <qmean/packing_statistics.hh>

namespace qmean {

class PackingPotential;
typedef boost::shared_ptr<PackingPotential> PackingPotentialPtr;
typedef MultiClassifier<float, int, int> PackingEnergies;


class DLLEXPORT_QMEAN PackingPotential : public PotentialBase, public impl::PackingPotentialImpl{
public:

  PackingPotential() { }
  
  static PackingPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& counts);

  static PackingPotentialPtr Create(PackingStatisticPtr stat, Real sigma, const String& reference_state);

  void Fill(PackingStatisticPtr stat, const String& reference_state);

  PotentialType GetType() const { return Packing; }

  void SetEnvironment(ost::mol::EntityView& env) { env_=env; }

  Real GetEnergy(atom::ChemType a, int count);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  impl::PackingOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;
    ds & energies_;
  }

private:
  
  PackingEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // PACKING_POTENTIAL_HH

