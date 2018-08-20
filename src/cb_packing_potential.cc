// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

