#include <qmean/interaction_potential.hh>


namespace qmean{

InteractionPotentialPtr InteractionPotential::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open interaction potential. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  InteractionPotentialPtr p(new InteractionPotential);
  ds >> *p;
  return p;
}

void InteractionPotential::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}

void InteractionPotential::OnInteraction(atom::ChemType a, atom::ChemType b, Real dist){
  energy_+=energies_.Get(a, b, dist);
  ++count_;
}

InteractionPotentialPtr InteractionPotential::Create(InteractionStatisticPtr s, Real sigma, const String& reference_state){

  InteractionPotentialPtr p(new InteractionPotential);
  p->opts_=s->GetOpts();
  p->opts_.sigma=sigma;
  p->energies_=InteractionEnergies(0.0, IntegralClassifier(atom::UNKNOWN, 0),
                                        IntegralClassifier(atom::UNKNOWN, 0),
                                        ContinuousClassifier(p->opts_.number_of_bins,
                                                             p->opts_.lower_cutoff,
                                                             p->opts_.upper_cutoff));

  p->Fill(s, reference_state);
  return p;
}

void InteractionPotential::Fill(InteractionStatisticPtr stat, const String& reference_state)
{
  typedef InteractionEnergies::IndexType Index;
  Real total_count = stat->GetTotalCount();

  boost::multi_array<Real,1> reference(boost::extents[opts_.number_of_bins]);

  if(reference_state=="classic"){

    for(int i=0;i<opts_.number_of_bins;++i){
      reference[i]=stat->GetCount(i)/total_count;
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in interaction potential!"
       << "implemented are: classic";
    throw io::IOException(ss.str());

  }
  
  for (int i=0; i<atom::UNKNOWN; ++i) {
    for (int j=0; j<atom::UNKNOWN; ++j) {
      Real sequence_count=stat->GetCount(atom::ChemType(i), atom::ChemType(j));
      for (int k=0; k<opts_.number_of_bins; ++k) {
        Real sequence_conformation_count=stat->GetCount(atom::ChemType(i), atom::ChemType(j),k);
        Real propensity=0.0;
        if(reference[k]>0 && sequence_count>0){
          propensity=sequence_conformation_count/sequence_count/reference[k];
        }
        Real e=log(1+opts_.sigma*sequence_count)-log(1+opts_.sigma*sequence_count*propensity);
        energies_.Set(Index(i, j, k), e);
      }
    }
  }
}

Real InteractionPotential::GetEnergy(atom::ChemType a, atom::ChemType b, Real dist){

  if(dist<opts_.lower_cutoff || dist>=opts_.upper_cutoff){
    return std::numeric_limits<Real>::quiet_NaN();
  }
  if(a==atom::UNKNOWN || b==atom::UNKNOWN){
    return std::numeric_limits<Real>::quiet_NaN();    
  }
  return energies_.Get(a, b, dist);
}

Real InteractionPotential::GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize){
  env_=env;
  return this->GetEnergy(target, normalize);
}

Real InteractionPotential::GetEnergy(mol::ResidueView& target_view, bool normalize){
  energy_=0.0;
  count_=0;
  target_view.Apply(*this);
  if(count_>0){
    if(normalize){
      return energy_/count_;
    }
    return energy_;
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

std::vector<Real> InteractionPotential::GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  env_=env;
  std::vector<Real> result;
  ost::mol::ResidueViewList res_list = target.GetResidueList();
  for(ost::mol::ResidueViewList::iterator it=res_list.begin(); it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it, normalize));
  }
  return result;
}

Real InteractionPotential::GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
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

}//namespace
