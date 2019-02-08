// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
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

#include <qmean/cbeta_impl.hh>

namespace qmean{ namespace impl{

void CBetaPotentialImpl::SetEnvironment(ost::mol::EntityView& env){

  env_.Clear();

  ost::mol::AtomView cb;
  ost::mol::AtomView ca;
  ost::mol::AtomView c;
  ost::mol::AtomView n;
  geom::Vec3 pos;
  int number;

  ost::mol::ResidueViewList res_list = env.GetResidueList();

  for(ost::mol::ResidueViewList::iterator it = res_list.begin(); it!=res_list.end(); ++it){

    ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(it->GetOneLetterCode());
    if(aa==ost::conop::XXX){
      continue;
    }

    if(!it->IsPeptideLinking()){
      continue;
    }

    cb = it->FindAtom("CB");
    number = it->GetNumber().GetNum();

    if(cb.IsValid()){
      pos = cb.GetPos();
      env_.Add(CBetaSpatialOrganizerItem(aa, pos, number, it->GetChain().GetName()), pos);
      continue;
    }

    ca = it->FindAtom("CA");
    n = it->FindAtom("N");
    c =it->FindAtom("C");

    if( !(ca.IsValid() && n.IsValid() && c.IsValid()) ){
      //cannot reconstruct cbeta!
      //TODO: make some noise
      continue;
    }

    pos = ost::mol::alg::CBetaPosition(n.GetPos(), ca.GetPos(), c.GetPos());
    env_.Add(CBetaSpatialOrganizerItem(aa, pos, number, it->GetChain().GetName()), pos);
  }
}

bool CBetaPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  if (!res.IsPeptideLinking()) {
    return false;
  }
  ost::conop::AminoAcid aa=ost::conop::OneLetterCodeToAminoAcid(res.GetOneLetterCode());
  if (aa==ost::conop::XXX) {
    return false;
  }

  ost::mol::AtomHandle cb = res.FindAtom("CB");
  int number = res.GetNumber().GetNum();
  geom::Vec3 pos;
  String chain_name = res.GetChain().GetName();

  if(cb.IsValid()){
    pos = cb.GetPos();
  }
  else{
    
    ost::mol::AtomHandle ca = res.FindAtom("CA");
    ost::mol::AtomHandle n = res.FindAtom("N");
    ost::mol::AtomHandle c =res.FindAtom("C");

    if( !(ca.IsValid() && n.IsValid() && c.IsValid()) ){
      //cannot reconstruct cbeta!
      //TODO: make some noise
      return false;
    }

    pos = ost::mol::alg::CBetaPosition(n.GetPos(), ca.GetPos(), c.GetPos());
  }

  std::vector<CBetaSpatialOrganizerItem> in_reach=env_.FindWithin(pos, opts_.upper_cutoff);

  for(std::vector<CBetaSpatialOrganizerItem>::iterator it=in_reach.begin();it!=in_reach.end();++it){
    if(std::abs(it->number-number)<opts_.sequence_sep){
      if(chain_name==it->chain_name){
        continue;
      }
    }
    Real dist=geom::Distance(it->pos, pos);
    if(dist<opts_.lower_cutoff || dist>=opts_.upper_cutoff){
      continue;
    }
    this->OnInteraction(aa, it->aa, dist);
  }
}

}}

