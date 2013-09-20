#ifndef INTERACTION_IMPL_HH
#define INTERACTION_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <qmean/atom_types.hh>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN InteractionOpts{
InteractionOpts() { }

InteractionOpts(Real l_cutoff, Real u_cutoff, int nob, int ssep, Real s=0.02):
             lower_cutoff(l_cutoff), upper_cutoff(u_cutoff), number_of_bins(nob),
             sequence_sep(ssep), sigma(s)
  { }

public:
  Real lower_cutoff;
  Real upper_cutoff;
  int number_of_bins;
  int sequence_sep;
  Real sigma;


  template <typename DS>
  void Serialize(DS& ds)
  {
    ds & lower_cutoff;
    ds & upper_cutoff;
    ds & number_of_bins;
    ds & sequence_sep;
    ds & sigma;
  }
};

class DLLEXPORT_QMEAN InteractionPotentialImpl : public ost::mol::EntityVisitor {
public:

  InteractionPotentialImpl() { }

  InteractionPotentialImpl(InteractionOpts& opts):
    opts_(opts) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual bool VisitAtom(const ost::mol::AtomHandle& at);

  virtual void OnInteraction(atom::ChemType a,
                             atom::ChemType b,
                             Real dist)=0;

protected:
  InteractionOpts opts_;
  ost::mol::EntityView  env_;

private:
  int actual_residue_number_;
  ost::conop::AminoAcid actual_aa_;
  String actual_chain_name_;

};

}} //namespace




#endif // MREDUCED_IMPL_HH