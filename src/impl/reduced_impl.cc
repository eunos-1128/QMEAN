// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
// Biozentrum - University of Basel
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <qmean/reduced_impl.hh>

namespace qmean { namespace impl {

void ReducedPotentialImpl::SetEnvironment(ost::mol::EntityView& env){

  env_.Clear();

  ost::mol::ResidueViewList res_list = env.GetResidueList();

  geom::Vec3 pos;
  geom::Vec3 axis;
  int number;
  String chain_name;

  for(ost::mol::ResidueViewList::iterator it=res_list.begin();it!=res_list.end();++it){

    ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(it->GetOneLetterCode());
    if(aa==ost::conop::XXX){
      continue;
    }

    if (!this->GetAxis(it->GetHandle(), axis)) {
      continue;
    }
    number = it->GetNumber().GetNum();
    chain_name = it->GetChain().GetName();
    ost::mol::AtomView ca = it->FindAtom("CA");
    pos = ca.GetPos();

    env_.Add(ReducedSpatialOrganizerItem(pos, axis, aa, number, chain_name), pos);
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

  geom::Vec3 axis;
  Real alpha;
  Real beta;
  Real gamma;
  Real dist;
  String chain_name = res.GetChain().GetName();

  int number=res.GetNumber().GetNum();
  if (!this->GetAxis(res, axis)) {
    return false;
  }

  //ca will be valid, otherwise we wouldn't have reached this point
  ost::mol::AtomHandle ca = res.FindAtom("CA");
  geom::Vec3 ca_pos = ca.GetPos();

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
    geom::Vec3 vec_to_other = it->pos - ca_pos;
    alpha = geom::Angle(axis, vec_to_other);
    beta = geom::Angle(it->axis, -vec_to_other);
    gamma = geom::DihedralAngle(ca_pos+axis,ca_pos,it->pos,it->pos+it->axis);
    this->OnInteraction(aa, it->aa, dist, alpha, beta, gamma);
  }
  return false;
}

bool ReducedPotentialImpl::GetAxis(const ost::mol::ResidueHandle& res,
                                   geom::Vec3& axis){
  ost::mol::AtomHandle ca=res.FindAtom("CA");
  ost::mol::AtomHandle n=res.FindAtom("N");
  ost::mol::AtomHandle c=res.FindAtom("C");
  if (!(ca.IsValid() && c.IsValid() && n.IsValid())) {
    return false;
  }

  geom::Vec3 n_ca = geom::Normalize(ca.GetPos()-n.GetPos());
  geom::Vec3 c_ca = geom::Normalize(ca.GetPos()-c.GetPos());

  axis = n_ca + c_ca;

  return true;
}

}}

