#ifndef REDUCED_IMPL_HH
#define REDUCED_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>
#include <ost/mol/spatial_organizer.hh>
#include <ost/mol/alg/construct_cbeta.hh>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN ReducedOpts{
ReducedOpts() { }

ReducedOpts(Real l_cutoff, Real u_cutoff, int nab, int ndab,
            int ndb, int ssep, Real s=0.02):
            lower_cutoff(l_cutoff), upper_cutoff(u_cutoff),
            num_angle_bins(nab),  num_dihedral_bins(ndab),
            num_dist_bins(ndb), sequence_sep(ssep), sigma(s) { 

  if(lower_cutoff >= upper_cutoff){
    throw std::runtime_error("Lower cutoff must be smaller than upper cutoff!");
  }

  if(lower_cutoff < 0.0){
    throw std::runtime_error("Lower cutoff must be larger or equal 0.0");
  }

  if(num_angle_bins <= 0){
    throw std::runtime_error("Number of angle bins must be larger 0!");
  }

  if(num_dihedral_bins <= 0){
    throw std::runtime_error("Number of dihedral bins must be larger 0!");
  }

  if(num_dist_bins <= 0){
    throw std::runtime_error("Number of distance bins must be larger 0!");
  }

  if(sequence_sep < 0){
    throw std::runtime_error("Sequence separation must be larger or equal 0!");
  }

  if(sigma <= 0.0){
    throw std::runtime_error("Sigma parameter must be larger 0.0!");
  }
}

public:

  Real lower_cutoff;
  Real upper_cutoff;
  int num_angle_bins;
  int num_dihedral_bins;
  int num_dist_bins;
  int sequence_sep;
  Real sigma;

  inline bool operator == (const ReducedOpts& opts) const {
    return lower_cutoff == opts.lower_cutoff && 
           upper_cutoff == opts.upper_cutoff &&
           num_angle_bins == opts.num_angle_bins && 
           num_dihedral_bins == opts.num_dihedral_bins &&
           num_dist_bins == opts.num_dist_bins && 
           sequence_sep == opts.sequence_sep &&
           sigma == opts.sigma;
  }

  inline bool operator != (const ReducedOpts& opts) const {
    return !(*this == opts);
  }



  template <typename DS>
  void Serialize(DS& ds)
  {
    ds & lower_cutoff;
    ds & upper_cutoff;
    ds & num_angle_bins;
    ds & num_dihedral_bins;
    ds & num_dist_bins;
    ds & sequence_sep;
    ds & sigma;
  }
};

struct DLLEXPORT_QMEAN ReducedSpatialOrganizerItem{

  ReducedSpatialOrganizerItem() { }

  ReducedSpatialOrganizerItem(geom::Vec3 p, geom::Vec3 ax, 
                              ost::conop::AminoAcid a, int n, const String& c_n):
                              pos(p), axis(ax), aa(a), number(n), chain_name(c_n) { }

  geom::Vec3 pos;
  geom::Vec3 axis;
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

  static bool GetAxis(const ost::mol::ResidueHandle& res,
                       geom::Vec3& axis);

  virtual void OnInteraction(ost::conop::AminoAcid aa_one,
                             ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta,
                             Real gamma)=0;

  void SetEnvironment(ost::mol::EntityView& env);

protected:
  ReducedOpts opts_;
  ost::mol::SpatialOrganizer<ReducedSpatialOrganizerItem, geom::Vec3>  env_;
};

}} //namespace




#endif // MREDUCED_IMPL_HH

