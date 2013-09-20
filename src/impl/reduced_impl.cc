#include <qmean/reduced_impl.hh>

namespace qmean { namespace impl {

void ReducedPotentialImpl::SetEnvironment(ost::mol::EntityView& env){

  env_.Clear();

  ost::mol::ResidueViewList res_list = env.GetResidueList();

  geom::Vec3 pos;
  geom::Vec3 cb_pos;
  geom::Vec3 axis;
  int number;

  for(ost::mol::ResidueViewList::iterator it=res_list.begin();it!=res_list.end();++it){

    ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(it->GetOneLetterCode());
    if(aa==ost::conop::XXX){
      continue;
    }

    if (!this->GetCAlphaCBetaPos(it->GetHandle(), pos, cb_pos)) {
      continue;
    }

    axis = cb_pos-pos;
    number = it->GetNumber().GetNum();
    env_.Add(ReducedSpatialOrganizerItem(pos, axis, aa, number, it->GetChain().GetName()), pos);
  }
}

bool ReducedPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){
  if (!res.IsPeptideLinking()) {
    return false;
  }
  ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(res.GetOneLetterCode());
  if (aa==ost::conop::XXX) {
    return false;
  }
  geom::Vec3 ca_pos;
  geom::Vec3 cb_pos;
  geom::Vec3 axis;
  Real angle;
  Real dist;
  String chain_name = res.GetChain().GetName();

  int number=res.GetNumber().GetNum();
  if (!this->GetCAlphaCBetaPos(res, ca_pos, cb_pos)) {
    return false;
  }

  axis = cb_pos-ca_pos;

  std::vector<ReducedSpatialOrganizerItem> within=env_.FindWithin(ca_pos,
                                             opts_.upper_cutoff);
  
  for(std::vector<ReducedSpatialOrganizerItem>::iterator it=within.begin();it!=within.end();++it){
    dist = geom::Distance(it->pos, ca_pos);
    if(dist<opts_.lower_cutoff || dist>=opts_.upper_cutoff){
      continue;
    }
    if(std::abs(it->number-number)<opts_.sequence_sep){
      if(chain_name==it->chain_name){
        continue;
      }
    }

    angle = fmod(geom::Angle(axis, it->cacb_axis),M_PI);

    if(angle==M_PI){
      angle-=0.0001;
    }

    this->OnInteraction(aa, it->aa, dist, angle);
  }
  return false;
}

bool ReducedPotentialImpl::GetCAlphaCBetaPos(const ost::mol::ResidueHandle& res,
                                             geom::Vec3& ca_pos,
                                             geom::Vec3& cb_pos){
  ost::mol::AtomHandle ca=res.FindAtom("CA");
  if (!ca.IsValid()) {
    return false;
  }
  ca_pos=ca.GetPos();
  ost::mol::AtomHandle cb=res.FindAtom("CB");
  if (cb.IsValid()) {
    cb_pos=cb.GetPos();
    return true;
  }
  ost::mol::AtomHandle n=res.FindAtom("N");
  ost::mol::AtomHandle c=res.FindAtom("C");
  if (!(ca.IsValid() && c.IsValid() && n.IsValid())) {
    return false;
  }
  cb_pos=ost::mol::alg::CBetaPosition(n.GetPos(), ca.GetPos(), c.GetPos());
  return true;
}

}}

