#ifndef INTERACTION_POTENTIAL_HH
#define INTERACTION_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/interaction_impl.hh>
#include <qmean/interaction_statistics.hh>

namespace qmean {

class InteractionPotential;
typedef boost::shared_ptr<InteractionPotential> InteractionPotentialPtr;
typedef MultiClassifier<float, int, int, float> InteractionEnergies;


class DLLEXPORT_QMEAN InteractionPotential : public PotentialBase, public impl::InteractionPotentialImpl{
public:

  InteractionPotential() { }
  
  static InteractionPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(atom::ChemType a, atom::ChemType b, Real dist);

  static InteractionPotentialPtr Create(InteractionStatisticPtr s, Real sigma, const String& reference_state);

  void Fill(InteractionStatisticPtr stat, const String& reference_state);

  PotentialType GetType() { return Interaction; }

  void SetEnvironment(ost::mol::EntityView& env) { env_=env; }

  Real GetEnergy(atom::ChemType a, atom::ChemType b, Real dist);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() { return count_; }

  impl::InteractionOpts& GetOpts() { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & energies_;
    ds & opts_;
  }

private:
  
  InteractionEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // INTERACTION_POTENTIAL_HH