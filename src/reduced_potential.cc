#include <qmean/reduced_potential.hh>

namespace qmean{

ReducedPotentialPtr ReducedPotential::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open reduced potential. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  ReducedPotentialPtr p(new ReducedPotential);
  ds >> *p;
  return p;
}

void ReducedPotential::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}

void ReducedPotential::OnInteraction(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                             Real dist, Real alpha, Real beta, Real gamma){

  energy_+=energies_.Get(aa_one, aa_two, dist, alpha, beta, gamma);
  ++count_;
}

ReducedPotentialPtr ReducedPotential::Create(ReducedStatisticPtr s, Real sigma, const String& reference_state){

  ReducedPotentialPtr p(new ReducedPotential);
  p->opts_=s->GetOpts();
  p->opts_.sigma=sigma;
  p->energies_=ReducedEnergies(0.0, IntegralClassifier(20, 0),
                                    IntegralClassifier(20, 0),
                                    ContinuousClassifier(p->opts_.num_dist_bins,
                                                         p->opts_.lower_cutoff,
                                                         p->opts_.upper_cutoff),
                                    ContinuousClassifier(p->opts_.num_angle_bins,
                                                         0.0, M_PI),
                                    ContinuousClassifier(p->opts_.num_angle_bins,
                                                         0.0, M_PI),
                                    ContinuousClassifier(p->opts_.num_dihedral_bins,
                                                         -M_PI, M_PI));

  p->Fill(s, reference_state);
  return p;
}

void ReducedPotential::Fill(ReducedStatisticPtr stats, const String& reference_state){

  typedef ReducedEnergies::IndexType Index;
  boost::multi_array<Real,4> reference(boost::extents[opts_.num_dist_bins][opts_.num_angle_bins][opts_.num_angle_bins][opts_.num_dihedral_bins]);

  Real total_count=stats->GetTotalCount();

  if(reference_state=="classic"){
    for(int i=0;i<opts_.num_dist_bins;++i){
      for(int j=0;j<opts_.num_angle_bins;++j){
        for(int k=0;k<opts_.num_angle_bins;++k){
          for(int l=0;l<opts_.num_dihedral_bins;++l){
            reference[i][j][k][l]=stats->GetCount(i,j,k,l)/total_count;
          }
        }
      }
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in reduced potential!"
       << "implemented are: classic";
    throw io::IOException(ss.str());
  }

  for (size_t i=0; i<ost::conop::XXX; ++i) {
    for (size_t j=0; j<ost::conop::XXX; ++j) {
      Real sequence_count = stats->GetCount(ost::conop::AminoAcid(i), ost::conop::AminoAcid(j));
      for (size_t k=0; k<opts_.num_dist_bins; ++k) {
        for (size_t l=0; l<opts_.num_angle_bins; ++l) {
          for (size_t m=0; m<opts_.num_angle_bins; ++m) {
            for (size_t n=0; n<opts_.num_dihedral_bins; ++n) {
              Real propensity = 0.0;
              Real sequence_conformation_count = stats->GetCount(ost::conop::AminoAcid(i), ost::conop::AminoAcid(j), k, l, m, n);
              if(sequence_count>0 && reference[k][l][m][n]>0){
                propensity=sequence_conformation_count/sequence_count/reference[k][l][m][n];
              }
              Real e=(log(1+opts_.sigma*sequence_count)-log(1+opts_.sigma*sequence_count*propensity));
              energies_.Set(Index(i, j, k, l, m, n), e);
            }
          }
        }
      }
    }
  }
}

Real ReducedPotential::GetEnergy(ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two,
                                    Real dist, Real alpha, Real beta, Real gamma){

  if(aa_one == ost::conop::XXX || aa_two == ost::conop::XXX){
    throw std::runtime_error("Cannot get energy from invalid amino acid!");
  }
  if(alpha >= M_PI || alpha < 0){
    throw std::runtime_error("alpha must be within [0,pi[");
  }
  if(beta >= M_PI || beta < 0){
    throw std::runtime_error("beta must be within [0,pi[");
  }
  if(gamma >= M_PI || gamma < -M_PI){
    throw std::runtime_error("gamma must be within [-pi,pi[");
  }

  if(dist<opts_.lower_cutoff || dist>=opts_.upper_cutoff){
    return std::numeric_limits<Real>::quiet_NaN();
  }

  return energies_.Get(aa_one, aa_two, dist, alpha, beta, gamma);
}

Real ReducedPotential::GetEnergy(ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize){
  this->SetEnvironment(env);
  return this->GetEnergy(target, normalize);
}

Real ReducedPotential::GetEnergy(mol::ResidueView& target_view, bool normalize){
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

std::vector<Real> ReducedPotential::GetEnergies(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  this->SetEnvironment(env);
  std::vector<Real> result;
  ost::mol::ResidueViewList res_list = target.GetResidueList();
  for(ost::mol::ResidueViewList::iterator it=res_list.begin(); it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it, normalize));
  }
  return result;
}

Real ReducedPotential::GetTotalEnergy(ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
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

}//namespace

