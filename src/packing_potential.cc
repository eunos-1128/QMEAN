#include <qmean/packing_potential.hh>

using namespace ost;
using namespace ost::conop;

namespace qmean {


PackingPotentialPtr PackingPotential::Load(const String& filename){
 
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  PackingPotentialPtr p(new PackingPotential);
  ds >> *p;
  return p;
}

void PackingPotential::OnSave(io::BinaryDataSink& ds){
  ds <<  *this;
}

void PackingPotential::OnInteraction(ost::conop::AminoAcid aa,
                                     const ost::mol::AtomHandleList& a_list,
                                     std::vector<int>& bins){

  ost::mol::AtomHandleList::const_iterator a_it = a_list.begin();
  std::vector<int>::iterator b_it = bins.begin();

  if(a_list.size()!=bins.size()){
    std::cout<<"sizes are not equal!"<<std::endl;
  }

  for(;a_it!=a_list.end();++a_it,++b_it){
    atom::ChemType a = qmean::GetAtomTypeByName(aa, a_it->GetName());
    if(*b_it!=*b_it){
      continue;
    }
    if(a!=atom::UNKNOWN){
      energy_+=energies_.Get(a, *b_it);
      ++count_;
    }
  }
}


PackingPotentialPtr PackingPotential::Create(PackingStatisticPtr stat, Real sigma, const String& reference_state){
  PackingPotentialPtr p(new PackingPotential);
  p->opts_=stat->GetOpts();
  p->opts_.sigma=sigma;
  p->energies_=PackingEnergies(0.0, IntegralClassifier(atom::UNKNOWN, 0),
                                    IntegralClassifier(int(floor(p->opts_.max_counts/p->opts_.bin_size))+1, 0));

  p->Fill(stat, reference_state);
  return p;
}

PackingPotentialPtr PackingPotential::Create(PackingStatisticPtr s1, PackingStatisticPtr s2, Real sigma, const String& reference_state, Real max_energy){

  PackingPotentialPtr p(new PackingPotential);

  if(s1->GetOpts() != s1->GetOpts()){
    throw io::IOException("The options for both statistics must be consistent!");
  }

  p->opts_=s1->GetOpts();
  p->opts_.sigma=sigma;
  p->energies_=PackingEnergies(0.0, IntegralClassifier(atom::UNKNOWN, 0),
                                    IntegralClassifier(int(floor(p->opts_.max_counts/p->opts_.bin_size))+1, 0));

  p->Fill(s1, s2, reference_state, max_energy);
  return p;
}

void PackingPotential::Fill(PackingStatisticPtr stat, const String& reference_state)
{
  
  typedef PackingEnergies::IndexType Index;
  Real total_count=stat->GetTotalCount();

  boost::multi_array<Real,1> reference(boost::extents[opts_.max_counts/opts_.bin_size+1]);

  if(reference_state=="classic"){
    for(int i=0;i<opts_.max_counts/opts_.bin_size+1;++i){
      reference[i]=stat->GetCount(i)/total_count;
    }
  }
  else if(reference_state=="uniform"){
    for(int i=0;i<opts_.max_counts/opts_.bin_size+1;++i){
      reference[i]=Real(1.0)/(floor(opts_.max_counts/opts_.bin_size)+1);
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in packing potential!"
       << "implemented are: classic, uniform";
    throw io::IOException(ss.str());

  }

  for (int i=0; i < atom::UNKNOWN; ++i) {
    Real sequence_count = stat->GetCount(atom::ChemType(i));

    for (int j=0; j<opts_.max_counts/opts_.bin_size+1; ++j) {
      Real sequence_conformation_count=stat->GetCount(atom::ChemType(i), j);
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

void PackingPotential::Fill(PackingStatisticPtr s1, PackingStatisticPtr s2, const String& reference_state, Real max_energy)
{
  
  typedef PackingEnergies::IndexType Index;

  Real total_count_one=s1->GetTotalCount();
  Real total_count_two=s2->GetTotalCount();

  boost::multi_array<Real,1> reference_one(boost::extents[opts_.max_counts/opts_.bin_size+1]);
  boost::multi_array<Real,1> reference_two(boost::extents[opts_.max_counts/opts_.bin_size+1]);

  if(reference_state=="classic"){
    for(int i=0;i<opts_.max_counts/opts_.bin_size+1;++i){
      reference_one[i]=s1->GetCount(i)/total_count_one;
      reference_two[i]=s2->GetCount(i)/total_count_two;
    }
  }
  else if(reference_state=="uniform"){
    for(int i=0;i<opts_.max_counts/opts_.bin_size+1;++i){
      reference_one[i]=Real(1.0)/(floor(opts_.max_counts/opts_.bin_size)+1);
      reference_two[i]=Real(1.0)/(floor(opts_.max_counts/opts_.bin_size)+1);
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in packing potential!"
       << "implemented are: classic, uniform";
    throw io::IOException(ss.str());

  }

  for (int i=0; i < atom::UNKNOWN; ++i) {
    Real sequence_count_one = s1->GetCount(atom::ChemType(i));
    Real sequence_count_two = s2->GetCount(atom::ChemType(i));

    for (int j=0; j<opts_.max_counts/opts_.bin_size+1; ++j) {
      Real sequence_conformation_count_one=s1->GetCount(atom::ChemType(i), j);
      Real sequence_conformation_count_two=s2->GetCount(atom::ChemType(i), j);

      Real propensity_sum=0.0;
      if(sequence_count_one>0 && reference_one[j]>0) {
        propensity_sum += sequence_conformation_count_one/sequence_count_one/reference_one[j];
      }
      if(sequence_count_two>0 && reference_two[j]>0) {
        propensity_sum += sequence_count_two*opts_.sigma*sequence_conformation_count_two/sequence_count_two/reference_two[j];
      }  
      Real e;
      if(propensity_sum == 0) e = max_energy;
      else e = std::min(max_energy,Real(log(1+opts_.sigma*sequence_count_two)-log(propensity_sum)));
      energies_.Set(Index(i,j),e);
    }
  }
}

Real PackingPotential::GetEnergy(atom::ChemType a, int count){
  return energies_.Get(a, this->GetBin(count));
}

Real PackingPotential::GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize){
  env_=env;
  return this->GetEnergy(target, normalize);
}

Real PackingPotential::GetEnergy(ost::mol::ResidueView& target, bool normalize){
  energy_=0;
  count_=0;
  target.Apply(*this);

  if(count_==0){
    return std::numeric_limits<Real>::quiet_NaN();
  }
  else{
    if(normalize){
      return energy_/count_;
    }
    return energy_;
  }
}

std::vector<Real> PackingPotential::GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  env_=env;
  std::vector<Real> result;
  ost::mol::ResidueViewList res_list = target.GetResidueList();
  for(ost::mol::ResidueViewList::iterator it=res_list.begin(); it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it, normalize));
  }
  return result;
}

Real PackingPotential::GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  energy_=0;
  count_=0;
  env_=env;
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

