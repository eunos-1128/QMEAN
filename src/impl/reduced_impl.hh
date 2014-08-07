#ifndef REDUCED_IMPL_HH
#define REDUCED_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <ost/mol/spatial_organizer.hh>
#include <ost/mol/alg/construct_cbeta.hh>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN ReducedOpts{
ReducedOpts(): lower_cutoff(0), upper_cutoff(0), num_angular_bins(0),
                num_dist_bins(0), sequence_sep(0), sigma(0)
  { }

ReducedOpts(Real l_cutoff, Real u_cutoff, int nab,
             int ndb, int ssep, Real s=0.02):
             lower_cutoff(l_cutoff), upper_cutoff(u_cutoff),
             num_angular_bins(nab), num_dist_bins(ndb),
             sequence_sep(ssep), sigma(s)
  { }

public:

  Real lower_cutoff;
  Real upper_cutoff;
  int num_angular_bins;
  int num_dist_bins;
  int sequence_sep;
  Real sigma;


  template <typename DS>
  void Serialize(DS& ds)
  {
    ds & lower_cutoff;
    ds & upper_cutoff;
    ds & num_angular_bins;
    ds & num_dist_bins;
    ds & sequence_sep;
    ds & sigma;
  }
};

struct DLLEXPORT_QMEAN ReducedSpatialOrganizerItem{

  ReducedSpatialOrganizerItem() { }

  ReducedSpatialOrganizerItem(geom::Vec3 p, geom::Vec3 ax, 
                       ost::conop::AminoAcid a, int n, const String& c_n):
                       pos(p), cacb_axis(ax), aa(a), number(n), chain_name(c_n) { }

  geom::Vec3 pos;
  geom::Vec3 cacb_axis;
  ost::conop::AminoAcid aa;
  int number;
  String chain_name;
};

class DLLEXPORT_QMEAN ReducedPotentialImpl : public ost::mol::EntityVisitor {
public:

  ReducedPotentialImpl(): env_(5) { }

  ReducedPotentialImpl(ReducedOpts& opts):
    opts_(opts), env_(5) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  static bool GetCAlphaCBetaPos(const ost::mol::ResidueHandle& res,
                                geom::Vec3& ca_pos,
                                geom::Vec3& cb_pos);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one,
                             ost::conop::AminoAcid aa_two,
                             Real dist, Real angle)=0;

  void SetEnvironment(ost::mol::EntityView& env);

protected:
  ReducedOpts opts_;
  ost::mol::SpatialOrganizer<ReducedSpatialOrganizerItem, geom::Vec3>  env_;
};

}} //namespace




#endif // MREDUCED_IMPL_HH

