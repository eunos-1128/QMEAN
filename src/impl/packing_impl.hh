#ifndef PACKING_IMPL_HH
#define PACKING_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <qmean/atom_types.hh>
#include <limits>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN PackingOpts {

    PackingOpts () { }

    PackingOpts(Real co, int mc, int bs): cutoff(co), max_counts(mc), bin_size(bs), sigma(0)
    { }
    PackingOpts(Real co, int mc, int bs, Real s): cutoff(co), max_counts(mc), bin_size(bs), sigma(s)
    { }

    Real cutoff;
    int max_counts;
    int bin_size;
    Real sigma;

    template <typename DS>
    void Serialize(DS& ds)
    {
      ds & cutoff;
      ds & max_counts;
      ds & bin_size;
      ds & sigma;
    }
};

class DLLEXPORT_QMEAN PackingPotentialImpl : public ost::mol::EntityVisitor {
public:

  PackingPotentialImpl() { }

  PackingPotentialImpl(PackingOpts& opts):
    opts_(opts) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             const ost::mol::AtomHandleList& a_list,
                             std::vector<int>& bins)=0;

protected:

  int GetBin(int count);

  PackingOpts opts_;
  ost::mol::EntityView  env_;
};

}} //namespace




#endif // PACKING_IMPL_HH