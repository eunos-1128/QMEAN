#ifndef TORSION_IMPL_HH
#define TORSION_IMPL_HH

#include <ost/mol/mol.hh>
#include <qmean/module_config.hh>
#include <ost/string_ref.hh>
#include <ost/conop/amino_acids.hh>

namespace qmean { namespace impl {

struct DLLEXPORT_QMEAN TorsionOpts{

public:
  TorsionOpts(): sigma(0.02) { }
  TorsionOpts(std::vector<String>& gi, std::vector<int>& nob, Real s=0.02);

  template <typename DS>
  void Serialize(DS& ds){

    if(ds.IsSource()){
      int num_of_groups;
      ds & num_of_groups;
      for(int i=0;i<6;++i){
        num_of_bins.push_back(0);
      }
      for(int i=0;i<num_of_groups;++i){
        group_identifier.push_back("");
      }
    }
    else{
      int num_of_groups=group_identifier.size();
      ds & num_of_groups; 
    }
    for(std::vector<int>::iterator it=num_of_bins.begin();it!=num_of_bins.end();++it){
      ds & *it;
    }
    for(std::vector<String>::iterator it=group_identifier.begin();it!=group_identifier.end();++it){
      ds & *it;
    }
  }
 
  std::vector<int> num_of_bins;
  std::vector<String> group_identifier;
  Real sigma;
};

class DLLEXPORT_QMEAN TorsionPotentialImpl : public ost::mol::EntityVisitor {

public:

  TorsionPotentialImpl()  { }

  TorsionPotentialImpl(TorsionOpts& opts): opts_(opts) { }

  virtual bool VisitResidue(const ost::mol::ResidueHandle& res);

  virtual void OnInteraction(const String& group_identifier,
                             std::vector<Real>& angls)=0;

  String FindStat(std::vector<String>& residues);

protected:
  TorsionOpts opts_;
};

}}//namespace


#endif // TORSION_IMPL_HH
