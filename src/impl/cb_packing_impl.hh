#ifndef CB_PACKING_IMPL_HH
#define CB_PACKING_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <ost/mol/spatial_organizer.hh>
#include <ost/mol/alg/construct_cbeta.hh>
#include <limits>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN CBPackingOpts {

  CBPackingOpts () { }

  CBPackingOpts(Real co, int mc, int bs): cutoff(co), max_counts(mc), bin_size(bs), sigma(0)
  { }
  CBPackingOpts(Real co, int mc, int bs, Real s): cutoff(co), max_counts(mc), bin_size(bs), sigma(s)
  { }



  inline bool operator == (const CBPackingOpts& opts) const {
    return max_counts == opts.max_counts && cutoff == opts.cutoff &&
           bin_size == opts.bin_size && sigma == opts.sigma;
  }

  inline bool operator != (const CBPackingOpts& opts) const {
    return !(*this == opts);
  }


  Real cutoff;
  int max_counts;
  int bin_size;
  Real sigma;
  template <typename DS>
  void Serialize(DS& ds){
    ds & cutoff;
    ds & max_counts;
    ds & bin_size;
    ds & sigma;
  }
};

struct DLLEXPORT_QMEAN CBPackingSpatialOrganizerItem{

  CBPackingSpatialOrganizerItem(ost::conop::AminoAcid aa, const String& chain_name, int num):
                       aa(aa), chain_name(chain_name), num(num) { }
  ost::conop::AminoAcid aa;
  String chain_name;
  int num;
};

class DLLEXPORT_QMEAN CBPackingPotentialImpl : public ost::mol::EntityVisitor {
public:

  CBPackingPotentialImpl(): env_(5) { }

  CBPackingPotentialImpl(CBPackingOpts& opts):
    opts_(opts), env_(5) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(ost::conop::AminoAcid aa,
                             int bin)=0;

  void SetEnvironment(ost::mol::EntityView& env);

protected:

  int GetBin(int count);

  CBPackingOpts opts_;
  ost::mol::SpatialOrganizer<CBPackingSpatialOrganizerItem, geom::Vec3>  env_;
};

}} //namespace




#endif // PACKING_IMPL_HH
