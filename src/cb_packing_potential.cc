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

#include <qmean/cb_packing_potential.hh>

using namespace ost;
using namespace ost::conop;

namespace qmean {


CBPackingPotentialPtr CBPackingPotential::Load(const String& filename){
 
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  CBPackingPotentialPtr p(new CBPackingPotential);
  ds >> *p;
  return p;
}

void CBPackingPotential::OnSave(io::BinaryDataSink& ds){
  ds <<  *this;
}

void CBPackingPotential::OnInteraction(ost::conop::AminoAcid aa,
                                       int bin){
  if(aa != ost::conop::XXX){
    energy_+=energies_.Get(CBPackingEnergies::IndexType(aa, bin));
    ++count_;
  } 
}

CBPackingPotentialPtr CBPackingPotential::Create(CBPackingStatisticPtr stat, Real sigma, const String& reference_state){
  CBPackingPotentialPtr p(new CBPackingPotential);
  p->opts_=stat->GetOpts();
  p->opts_.sigma=sigma;
  p->energies_=CBPackingEnergies(0.0, IntegralClassifier(ost::conop::XXX, 0),
                                      IntegralClassifier(p->opts_.max_counts+1, 0));

  p->Fill(stat, reference_state);
  return p;
}

void CBPackingPotential::Fill(CBPackingStatisticPtr stat, const String& reference_state)
{
  
  typedef CBPackingEnergies::IndexType Index;
  Real total_count=stat->GetTotalCount();

  boost::multi_array<Real,1> reference(boost::extents[opts_.max_counts+1]);

  if(reference_state=="classic"){
    for(int i=0;i<opts_.max_counts+1;++i){
      reference[i]=stat->GetCount(i)/total_count;
    }
  }
  else if(reference_state=="uniform"){
    for(int i=0;i<opts_.max_counts+1;++i){
      reference[i]=Real(1.0)/(opts_.max_counts+1);
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in cb packing potential!"
       << "implemented are: classic, uniform";
    throw io::IOException(ss.str());
  }

  for (int i=0; i < ost::conop::XXX; ++i) {
    Real sequence_count = stat->GetCount(ost::conop::AminoAcid(i));

    for (int j=0; j<opts_.max_counts+1; ++j) {
      Real sequence_conformation_count=stat->GetCount(ost::conop::AminoAcid(i), j);
      Real propensity=0.0;
      if (sequence_count>0 && reference[j]>0) {
        propensity=sequence_conformation_count/sequence_count/reference[j];
      }
            
      Real e=log(1+opts_.sigma*sequence_count)-
               log(1+opts_.sigma*sequence_count*propensity);
      energies_.Set(Index(i,j),e);
    }
  }
}

Real CBPackingPotential::GetEnergy(ost::conop::AminoAcid aa, int count){

  if(aa == ost::conop::XXX){
    throw std::runtime_error("Cannot get energy for invalid amino acid!");
  }

  return energies_.Get(CBPackingEnergies::IndexType(aa, this->GetBin(count)));
}

Real CBPackingPotential::GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env){
  this->SetEnvironment(env);
  return this->GetEnergy(target);
}

Real CBPackingPotential::GetEnergy(ost::mol::ResidueView& target){
  energy_=0;
  count_=0;
  target.Apply(*this);

  if(count_==0){
    return std::numeric_limits<Real>::quiet_NaN();
  }
  else{
    return energy_;
  }
}

std::vector<Real> CBPackingPotential::GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env){
  this->SetEnvironment(env);
  std::vector<Real> result;
  ost::mol::ResidueViewList res_list = target.GetResidueList();
  for(ost::mol::ResidueViewList::iterator it=res_list.begin(); it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it));
  }
  return result;
}

Real CBPackingPotential::GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  energy_=0;
  count_=0;
  this->SetEnvironment(env);
  target.Apply(*this);
  if(count_>0){
    if(normalize){
      return energy_/count_;
    }
    return energy_;
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

}

