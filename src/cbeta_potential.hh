#ifndef CBETA_POTENTIAL_HH
#define CBETA_POTENTIAL_HH

#include <qmean/potential_base.hh>
#include <qmean/cbeta_impl.hh>
#include <qmean/cbeta_statistics.hh>


namespace qmean {

class CBetaPotential;
typedef boost::shared_ptr<CBetaPotential> CBetaPotentialPtr;
typedef MultiClassifier<float, int, int, float> CBetaEnergies;


class DLLEXPORT_QMEAN CBetaPotential : public PotentialBase, public impl::CBetaPotentialImpl{
public:

  CBetaPotential() { }
  
  static CBetaPotentialPtr Load(const String& file_name);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist);

  static CBetaPotentialPtr Create(CBetaStatisticPtr s, Real sigma, const String& reference_state);

  static CBetaPotentialPtr Create(CBetaStatisticPtr s1, CBetaStatisticPtr s2, Real sigma, const String& reference_state, Real max_energy);

  void Fill(CBetaStatisticPtr stat, const String& reference_state);

  void Fill(CBetaStatisticPtr s1, CBetaStatisticPtr s2, const String& reference_state, Real max_energy);

  PotentialType GetType() { return CBeta; }

  Real GetEnergy(ost::conop::AminoAcid a, ost::conop::AminoAcid, Real dist);

  Real GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize);

  Real GetEnergy(ost::mol::ResidueView& target, bool normalize);

  std::vector<Real> GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  Real GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize);

  uint32_t GetCounts() { return count_; }

  impl::CBetaOpts GetOpts() { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & energies_;
    ds & opts_;
  }

private:


  
  CBetaEnergies energies_;
  Real energy_;
  uint32_t count_;
};

}//namespace
#endif // CBETA_POTENTIAL_HH

