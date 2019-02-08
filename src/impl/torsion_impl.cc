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

#include <qmean/torsion_impl.hh>

namespace qmean { namespace impl{ 


TorsionOpts::TorsionOpts(std::vector<String>& gi, std::vector<int>& nob, Real s): sigma(s){

  if(nob.size()!=6){
    throw std::runtime_error("number of bin list must contain 6 values!");
  } 

  for(int i = 0; i < 6; ++i){
    if(nob[i] <= 0){
      throw std::runtime_error("All torsion angles must have at least one bin!");
    }
  }

  for(std::vector<String>::iterator it=gi.begin();it!=gi.end();++it){
    ost::StringRef gid(it->c_str(), it->length());
    if(gid.split('-').size()!=3){
      std::stringstream ss;
      ss<<"GroupIDs must have the form 'AANames-AANames-AANames' ";
      ss<<"whereas AANames are three letter abbreviations of amino acids ";
      ss<<"separated by ','. The keyword 'all' matches all amino acids. ";
      ss<<"e.g. 'all-VAL,LEU,ILE-all";
      throw std::runtime_error(ss.str());      
    }
  }

  num_of_bins=nob;
  group_identifier=gi;

  std::vector<String> single_ids(3);
  String current_identifier;
  for(int i = 0; i < ost::conop::XXX; ++i){
    single_ids[0] = ost::conop::AminoAcidToResidueName(ost::conop::AminoAcid(i));

    for(int j = 0; j < ost::conop::XXX; ++j){
      single_ids[1] = ost::conop::AminoAcidToResidueName(ost::conop::AminoAcid(j));

      for(int k = 0; k < ost::conop::XXX; ++k){
        single_ids[2] = ost::conop::AminoAcidToResidueName(ost::conop::AminoAcid(k));
        current_identifier = this->FindStat(single_ids);
        if(current_identifier == ""){
          std::stringstream ss;
          ss << "Amino acid triplet \""<<single_ids[0]<<", "<<single_ids[1]<<", ";
          ss << single_ids[2]<<"\" is not covered by any of the provided group ids!";
          throw std::runtime_error(ss.str());
        }
      }
    }
  }     
}

String TorsionOpts::FindStat(const std::vector<String>& residues){
  bool match;

  for(std::vector<String>::iterator it=group_identifier.begin();it!=group_identifier.end();++it){

    match=true;
    ost::StringRef gid(it->c_str(), it->length());
    std::vector<ost::StringRef> single_ids=gid.split('-'); //split ids

    //iterate over all three positions
    for(int i=0;i<3;++i){
      //"all" matches every residue
      if(single_ids[i].str().find("all")!=(-1)){
        continue;
      }
      //check, whether currrent residue matches current id position
      if(single_ids[i].str().find(residues[i])==-1){
        match=false;
        break;
      }
    }

    if(match){
      //FIRST match gets returned
      return *it;
    }
  }
  String empty("");
  return empty;
}

bool TorsionPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  std::vector<Real> angles;
  std::vector<String> r_names;

  if(res.GetPrev().IsValid()){
    ost::mol::ResidueHandle prev=res.GetPrev();
    if(!prev.IsPeptideLinking()){
      return false;
    }
    if(prev.GetOneLetterCode()==ost::conop::XXX){
      return false;
    }
    if(!prev.GetPhiTorsion().IsValid() || !prev.GetPsiTorsion().IsValid()){
      return false;
    }
    angles.push_back(std::min(prev.GetPhiTorsion().GetAngle(),Real(M_PI-0.0001)));
    angles.push_back(std::min(prev.GetPsiTorsion().GetAngle(),Real(M_PI-0.0001)));
    r_names.push_back(prev.GetName());
  }
  else{
    return false;
  }

  if(res.IsValid()){
    if(!res.IsPeptideLinking()){
      return false;
    }
    if(res.GetOneLetterCode()==ost::conop::XXX){
      return false;
    }
    if(!res.GetPhiTorsion().IsValid() || !res.GetPsiTorsion().IsValid()){
      return false;
    }
    angles.push_back(std::min(res.GetPhiTorsion().GetAngle(),Real(M_PI-0.0001)));
    angles.push_back(std::min(res.GetPsiTorsion().GetAngle(),Real(M_PI-0.0001)));
    r_names.push_back(res.GetName());
  }
  else{
    return false;
  }

  if(res.GetNext().IsValid()){
    ost::mol::ResidueHandle next=res.GetNext();
    if(!next.IsPeptideLinking()){
      return false;
    }
    if(next.GetOneLetterCode()==ost::conop::XXX){
      return false;
    }
    if(!next.GetPhiTorsion().IsValid() || !next.GetPsiTorsion().IsValid()){
      return false;
    }
    angles.push_back(std::min(next.GetPhiTorsion().GetAngle(),Real(M_PI-0.0001)));
    angles.push_back(std::min(next.GetPsiTorsion().GetAngle(),Real(M_PI-0.0001)));
    r_names.push_back(next.GetName());
  }
  else{
    return false;
  }

  String s = FindStat(r_names);

  if(s==""){
    return false;
  }

  this->OnInteraction(s,angles);

  return false;
}

String TorsionPotentialImpl::FindStat(std::vector<String>& residues){
  return opts_.FindStat(residues);
}

}}//namespaces

