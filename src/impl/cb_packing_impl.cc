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

#include <qmean/cb_packing_impl.hh>

namespace qmean{ namespace impl{

bool CBPackingPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  ost::conop::AminoAcid aa=ost::conop::ResidueToAminoAcid(res);
  if (aa==ost::conop::XXX){
    return false;
  }

  geom::Vec3 pos;
  ost::mol::AtomHandle cb;

  cb = res.FindAtom("CB");

  if(cb.IsValid()){
    pos = cb.GetPos();
  }
  else{
    ost::mol::AtomHandle ca = res.FindAtom("CA");
    ost::mol::AtomHandle c = res.FindAtom("C");
    ost::mol::AtomHandle n = res.FindAtom("N");
    if( !(ca.IsValid() && n.IsValid() && c.IsValid()) ){
      //cannot reconstruct cbeta!
      //TODO: make some noise
      return false;
    }
    pos = ost::mol::alg::CBetaPosition(n.GetPos(), ca.GetPos(), c.GetPos());
  }

  std::vector<CBPackingSpatialOrganizerItem> in_reach = env_.FindWithin(pos, opts_.cutoff);

  int count = 0;
  String chain_name = res.GetChain().GetName();
  int num = res.GetNumber().GetNum();

  for(std::vector<CBPackingSpatialOrganizerItem>::iterator i = in_reach.begin();
      i != in_reach.end(); ++i){
    if(i->chain_name == chain_name && i->num == num) continue;
    ++count;
  }

  this->OnInteraction(aa, this->GetBin(count));

  return false;
}

void CBPackingPotentialImpl::SetEnvironment(ost::mol::EntityView& env){

  env_.Clear();

  ost::mol::AtomView cb;
  ost::mol::AtomView ca;
  ost::mol::AtomView c;
  ost::mol::AtomView n;
  geom::Vec3 pos;
  String chain_name;
  int residue_number;

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
    chain_name = it->GetChain().GetName();
    residue_number = it->GetNumber().GetNum();

    if(cb.IsValid()){
      pos = cb.GetPos();
      env_.Add(CBPackingSpatialOrganizerItem(aa, chain_name, residue_number), pos);
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
    env_.Add(CBPackingSpatialOrganizerItem(aa, chain_name, residue_number), pos);
  }
}

int CBPackingPotentialImpl::GetBin(int count){
  int bin = std::min(count, opts_.max_counts);
  return bin;
}

}}

