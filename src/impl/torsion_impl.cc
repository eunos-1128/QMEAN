#include <qmean/torsion_impl.hh>

namespace qmean { namespace impl{ 


TorsionOpts::TorsionOpts(std::vector<String>& gi, std::vector<int>& nob, Real s): sigma(s){

  if(nob.size()!=6){
    throw std::runtime_error("number of bin list must contain 6 values!");
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

  this->OnInteraction(FindStat(r_names),angles);

  return false;
}

String TorsionPotentialImpl::FindStat(std::vector<String>& residues){

  bool match;

  for(std::vector<String>::iterator it=opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){

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

}}//namespaces

