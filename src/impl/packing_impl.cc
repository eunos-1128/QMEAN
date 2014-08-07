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
  int bin=int(std::floor((std::min(count, opts_.max_counts))/opts_.bin_size));
  return bin;
}

}}

