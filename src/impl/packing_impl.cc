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

#include <qmean/packing_impl.hh>

namespace qmean{ namespace impl{

bool PackingPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  ost::conop::AminoAcid aa=ost::conop::ResidueToAminoAcid(res);
  if (aa==ost::conop::XXX){
    return false;
  }

  ost::mol::AtomHandleList a_list=res.GetAtomList();
  std::vector<int> bins;
  int count;

  for(ost::mol::AtomHandleList::const_iterator it=a_list.begin();it!=a_list.end();++it){

    if(qmean::GetAtomTypeByName(aa,it->GetName())==atom::UNKNOWN){
      bins.push_back(std::numeric_limits<int>::quiet_NaN());
      continue;
    }

    ost::mol::AtomViewList in_reach=env_.FindWithin(it->GetPos(), opts_.cutoff);
    count=0; 

    for(ost::mol::AtomViewList::iterator ite=in_reach.begin(); ite!=in_reach.end(); ++ite){
      if(ite->GetHandle().GetResidue().GetName()=="MEM"){
        ++count;
        continue;
      }
      if(ost::conop::ResidueToAminoAcid(ite->GetHandle().GetResidue())==ost::conop::XXX){
         continue;
      }
      if(qmean::GetAtomTypeByName(ost::conop::ResidueToAminoAcid(ite->GetHandle().GetResidue()),ite->GetName())==atom::UNKNOWN){
        continue;
      }
      if(ite->GetHandle().GetResidue()!=res){
        ++count; 
      } 
    }
    bins.push_back(this->GetBin(count));
  }


  this->OnInteraction(aa, a_list, bins);

  return false;

}

int PackingPotentialImpl::GetBin(int count){
  int bin = std::min(count, opts_.max_counts);
  return bin;
}

}}

