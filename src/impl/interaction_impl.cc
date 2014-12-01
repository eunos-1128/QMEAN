#include <qmean/interaction_impl.hh>

namespace qmean{ namespace impl{

bool InteractionPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  if (!res.IsPeptideLinking()) {
    return false;
  }
  ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(res.GetOneLetterCode());
  if (aa==ost::conop::XXX) {
    return false;
  }

  actual_residue_number_=res.GetNumber().GetNum();
  actual_aa_=aa;
  actual_chain_name_=res.GetChain().GetName();

  return true;
}

bool InteractionPotentialImpl::VisitAtom(const ost::mol::AtomHandle& at){

  atom::ChemType type_a=qmean::GetAtomTypeByName(actual_aa_, at.GetName()); 

  if(type_a==atom::UNKNOWN){
    return false;
  }

  atom::ChemType type_b; 
  int b_residue_number;
  ost::conop::AminoAcid b_aa;

  geom::Vec3 pos = at.GetPos();
  Real dist;

  ost::mol::AtomViewList in_reach = env_.FindWithin(pos, opts_.upper_cutoff);

  for(ost::mol::AtomViewList::iterator it=in_reach.begin();it!=in_reach.end();++it){

    b_aa=ost::conop::OneLetterCodeToAminoAcid(it->GetResidue().GetOneLetterCode());
    if(b_aa==ost::conop::XXX){
      continue;
    }

    b_residue_number=it->GetResidue().GetNumber().GetNum();
    if(abs(b_residue_number-actual_residue_number_)<opts_.sequence_sep){
      if(it->GetResidue().GetChain().GetName()==actual_chain_name_){
        continue;
      }
    }

    type_b=qmean::GetAtomTypeByName(b_aa, it->GetName());
    if(type_b==atom::UNKNOWN){
      continue;
    }

    dist=geom::Distance(pos, it->GetPos());

    if(dist<opts_.lower_cutoff || dist>=opts_.upper_cutoff){
      continue;
    }

    this->OnInteraction(type_a, type_b, dist);

  }
}

}} //namespaces

