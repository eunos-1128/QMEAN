#ifndef HBOND_IMPL_HH
#define HBOND_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/mol/spatial_organizer.hh>
#include <limits>
#include <vector>

namespace qmean { namespace impl {

  
struct DLLEXPORT_QMEAN HBondOpts{
public:
  HBondOpts () { }

  HBondOpts(Real d_mi, Real d_ma, Real alpha_mi, Real alpha_ma,
            Real beta_mi, Real beta_ma, Real gamma_mi, Real gamma_ma, 
            int d_b, int alpha_b, int beta_b, 
            int gamma_b, int ss): d_min(d_mi), d_max(d_ma), 
                                  alpha_min(alpha_mi),alpha_max(alpha_ma),
                                  beta_min(beta_mi),beta_max(beta_ma),
                                  gamma_min(gamma_mi),gamma_max(gamma_ma),
                                  d_bins(d_b),alpha_bins(alpha_b),beta_bins(beta_b),
                                  gamma_bins(gamma_b), seq_sep(ss) { }

  inline bool operator == (const HBondOpts& opts) const {
    return d_min == opts.d_min && d_max == opts.d_max &&
           alpha_min == opts.alpha_min && alpha_max == opts.alpha_max &&
           beta_min == opts.beta_min && beta_max == opts.beta_max &&
           gamma_min == opts.gamma_min && gamma_max == opts.gamma_max &&
           d_bins == opts.d_bins && alpha_bins == opts.alpha_bins &&
           beta_bins == opts.beta_bins && gamma_bins == opts.gamma_bins &&
           seq_sep == opts.seq_sep;
  }

  inline bool operator != (const HBondOpts& opts) const {
    return !(*this == opts);
  }

  Real d_min;
  Real d_max;
  Real alpha_min;
  Real alpha_max;
  Real beta_min;
  Real beta_max;
  Real gamma_min;
  Real gamma_max;
  int d_bins;
  int alpha_bins;
  int beta_bins;
  int gamma_bins;
  int seq_sep;
  std::vector<Real> abs_expectation_values;

  template <typename DS>
  void Serialize(DS& ds){
    ds & d_min;
    ds & d_max;
    ds & alpha_min;
    ds & alpha_max;
    ds & beta_min;
    ds & beta_max;
    ds & gamma_min;
    ds & gamma_max;
    ds & d_bins;
    ds & alpha_bins;
    ds & beta_bins;
    ds & gamma_bins;
    ds & seq_sep;
    if(ds.IsSource()){
      int num_elements;
      ds & num_elements;
      for(int i = 0; i < num_elements; ++i){
        abs_expectation_values.push_back(0);
        ds & abs_expectation_values[i];
      }
    }
    else{
      int num_elements = abs_expectation_values.size();
      ds & num_elements;
      for(uint i = 0; i < abs_expectation_values.size(); ++i){
        ds & abs_expectation_values[i];
      }
    }
  }
};

struct DLLEXPORT_QMEAN HBondSpatialOrganizerItem{

  HBondSpatialOrganizerItem(geom::Vec3& n, geom::Vec3& ca, geom::Vec3& c,
                            geom::Vec3& o, geom::Vec3& h, int s,
                            String c_n, int nu):n_pos(n),
                                                ca_pos(ca),
                                                c_pos(c),
                                                o_pos(o),
                                                h_pos(h),
                                                state(s), 
                                                chain_name(c_n), 
                                                num(nu) { }
  geom::Vec3 n_pos;
  geom::Vec3 ca_pos;
  geom::Vec3 c_pos;
  geom::Vec3 o_pos;
  geom::Vec3 h_pos;
  int state;
  String chain_name;
  int num;
};

class DLLEXPORT_QMEAN HBondPotentialImpl : public ost::mol::EntityVisitor {
public:

  HBondPotentialImpl(): env_(5) { }

  HBondPotentialImpl(HBondOpts& opts):
    opts_(opts), env_(5) { }


  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(int state,Real d,Real alpha,
                             Real beta,Real gamma)=0;

  void SetEnvironment(const ost::mol::EntityView& env);

protected:

  HBondOpts opts_;
  ost::mol::SpatialOrganizer<HBondSpatialOrganizerItem, geom::Vec3>  env_;
  bool unique_bonds_;

private:
  int GetState(const ost::mol::ResidueHandle& res);

  void ConstructHPos(Real bond_length, Real angle, Real dihedral, 
                     geom::Vec3& A, geom::Vec3& B, geom::Vec3& C,
                     geom::Vec3& pos);
};

}} //namespace

#endif // HBOND_IMPL_HH
