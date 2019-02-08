#include <qmean/hbond_potential.hh>

namespace qmean {

HBondPotentialPtr HBondPotential::Load(const String& filename){
 
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open potential. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  HBondPotentialPtr p(new HBondPotential);
  ds >> *p;
  return p;
}

void HBondPotential::OnSave(io::BinaryDataSink& ds){
  ds <<  *this;
}

void HBondPotential::OnInteraction(int state, Real d, Real alpha, Real beta, Real gamma){
  energy_ += energies_.Get(state,d,alpha,beta,gamma);
}


HBondPotentialPtr HBondPotential::Create(HBondStatisticPtr stat){
  HBondPotentialPtr p(new HBondPotential);
  p->opts_=stat->GetOpts();
  p->energies_ = HBondEnergies(0.0, IntegralClassifier(3, 0),
                                    ContinuousClassifier(p->opts_.d_bins,p->opts_.d_min,p->opts_.d_max),
                                    ContinuousClassifier(p->opts_.alpha_bins,p->opts_.alpha_min,p->opts_.alpha_max),
                                    ContinuousClassifier(p->opts_.beta_bins,p->opts_.beta_min,p->opts_.beta_max),
                                    ContinuousClassifier(p->opts_.gamma_bins,p->opts_.gamma_min,p->opts_.gamma_max));

  p->Fill(stat);
  return p;
}

void HBondPotential::Fill(HBondStatisticPtr stat){
  
  typedef HBondEnergies::IndexType Index;

  double ref_value = 1.0 / (opts_.d_bins * opts_.alpha_bins * opts_.beta_bins * opts_.gamma_bins);
  double state_count;
  double count;
  double frac;
  Real e;
  double expectation_value;

  for (int i=0; i < 3; ++i) {
    state_count = stat->GetCount(i);
    if(state_count == 0.0){
      std::stringstream ss;
      ss <<"Cannot create potential function from state " << i;
      ss <<" without counts in histogram!";
      throw std::runtime_error(ss.str());
    }
    expectation_value = 0.0;
    for(uint j = 0; j < opts_.d_bins; ++j){
      for(uint k = 0; k < opts_.alpha_bins; ++k){
        for(uint l = 0; l < opts_.beta_bins; ++l){
          for(uint m = 0; m < opts_.gamma_bins; ++m){
            count = stat->GetCount(i,j,k,l,m);
            frac = std::max((count/state_count)/ref_value,1.0);
            e = -log(frac);
            energies_.Set(Index(i,j,k,l,m),e);
            expectation_value += count/state_count*e;
          }
        }
      }
    }
    expectation_value = std::abs(expectation_value);
    for(uint j = 0; j < opts_.d_bins; ++j){
      for(uint k = 0; k < opts_.alpha_bins; ++k){
        for(uint l = 0; l < opts_.beta_bins; ++l){
          for(uint m = 0; m < opts_.gamma_bins; ++m){
            e = energies_.Get(Index(i,j,k,l,m));
            energies_.Set(Index(i,j,k,l,m),e/expectation_value);            
          }
        }
      }
    }
  }
}

Real HBondPotential::GetEnergy(int state, float d, float alpha, float beta, float gamma){
  if(d<opts_.d_min || d>=opts_.d_max) return 0.0;
  if(alpha<opts_.alpha_min || alpha>=opts_.alpha_max) return 0.0;
  if(beta<opts_.beta_min || beta>=opts_.beta_max) return 0.0;
  if(gamma<opts_.gamma_min || gamma>=opts_.gamma_max) return 0.0;
  if(state < 0 || state >= 3) return 0.0;
  return energies_.Get(state, d, alpha, beta, gamma);
}

Real HBondPotential::GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env){
  this->SetEnvironment(env);
  return this->GetEnergy(target);
}

Real HBondPotential::GetEnergy(ost::mol::ResidueView& target){
  energy_=0;
  unique_bonds_ = false;
  target.Apply(*this);
  return energy_;
}

std::vector<Real> HBondPotential::GetEnergies(ost::mol::EntityView& view){
  this->SetEnvironment(view);
  std::vector<Real> result;
  ost::mol::ResidueViewList res_list = view.GetResidueList();
  for(ost::mol::ResidueViewList::iterator it=res_list.begin(); it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it));
  }
  return result;
}

Real HBondPotential::GetTotalEnergy(ost::mol::EntityView& view, bool normalize){
  energy_=0;
  int count = view.GetResidueCount();
  this->SetEnvironment(view);
  unique_bonds_ = true;
  view.Apply(*this);
  if(count > 0){
    if(normalize){
      return energy_/count;
    }
    return energy_;
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

} //namespace
