#include <qmean/torsion_potential.hh>

using namespace ost;

namespace qmean{

TorsionPotentialPtr TorsionPotential::Create(TorsionStatisticPtr& stat, Real sigma, const String& reference_state){

  TorsionPotentialPtr p(new TorsionPotential);

  p->opts_=stat->GetOpts();
  p->opts_.sigma=sigma;

  for(int i=0;i<p->opts_.group_identifier.size();++i){
    p->energies_[p->opts_.group_identifier[i]]=TorsionEnergies(0.0,
                                                  ContinuousClassifier(p->opts_.num_of_bins[0], -M_PI, M_PI),
                                                  ContinuousClassifier(p->opts_.num_of_bins[1], -M_PI, M_PI),
                                                  ContinuousClassifier(p->opts_.num_of_bins[2], -M_PI, M_PI),
                                                  ContinuousClassifier(p->opts_.num_of_bins[3], -M_PI, M_PI),
                                                  ContinuousClassifier(p->opts_.num_of_bins[4], -M_PI, M_PI),
                                                  ContinuousClassifier(p->opts_.num_of_bins[5], -M_PI, M_PI));
  }

  p->Fill(stat, reference_state);

  return p;
}

TorsionPotentialPtr TorsionPotential::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open torsion potential. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  TorsionPotentialPtr p(new TorsionPotential);
  ds >> *p;
  return p;
}

void TorsionPotential::OnSave(io::BinaryDataSink& ds){
  ds << *this;
}


void TorsionPotential::OnInteraction(const String& group_identifier,
                                        std::vector<Real>& angles){

  energy_+=energies_[group_identifier].Get(angles[0],angles[1],angles[2],
                                           angles[3],angles[4],angles[5]);
  ++count_;
}

std::vector<Real> TorsionPotential::GetEnergies(ost::mol::EntityView& target){

  std::vector<Real> result;
  ost::mol::ResidueViewList res_list=target.GetResidueList();

  for(ost::mol::ResidueViewList::iterator it=res_list.begin();it!=res_list.end();++it){
    result.push_back(this->GetEnergy(*it));
  }

  return result;
}

Real TorsionPotential::GetEnergy(ost::mol::ResidueView& target_view){
  energy_=0;
  count_=0;
  target_view.Apply(*this);
  if(count_==1){
    return energy_;
  }
  else{
    return std::numeric_limits<Real>::quiet_NaN();
  }
}

Real TorsionPotential::GetEnergy(const String& group_id, std::vector<Real>& angles){

  if(angles.size()!=6){
    std::stringstream ss;
    ss<<"requires 6 angles, ";
    ss<<angles.size();
    ss<<" provided.";
    throw std::runtime_error(ss.str());
  }

  energy_=0;
  count_=0;

  std::vector<Real> new_angles(angles);
  Real new_angle;

  if(energies_.find(group_id) == energies_.end()){
    return std::numeric_limits<Real>::quiet_NaN();
  }

  //let's make sure the angles are in a valid range
  for(uint i = 0; i < 6; ++i){
    while (new_angles[i] >= M_PI) new_angles[i] -= 2*M_PI;
    while (new_angles[i] < -M_PI) new_angles[i] += 2*M_PI;
  }

  this->OnInteraction(group_id, new_angles);

  if(count_==1){
    return energy_;
  }
  else{
    std::cout<<"no energy_counts"<<std::endl;
    return std::numeric_limits<Real>::quiet_NaN();
  }


}

Real TorsionPotential::GetEnergy(std::vector<String>& residue_names, std::vector<Real>& angles){

  if(residue_names.size()!=3){
    std::stringstream ss;
    ss<<"requires 3 residue names, only ";
    ss<<residue_names.size();
    ss<<" provided.";
    throw std::runtime_error(ss.str());
  }

  String group_id=this->FindStat(residue_names);

  return this->GetEnergy(group_id,angles);
}

Real TorsionPotential::GetTotalEnergy(ost::mol::EntityView& target, bool normalize){
  energy_=0;
  count_=0;
  target.Apply(*this);
  if(count_>0){
    if(normalize){
      return energy_/count_;;
    }
    return energy_;
  }
  return std::numeric_limits<Real>::quiet_NaN();
}

void TorsionPotential::Fill(TorsionStatisticPtr stat, const String& reference_state){

  typedef TorsionEnergies::IndexType Index;
  Real total_count=stat->GetTotalCount();

  if(total_count==0){
    LOG_ERROR("There are no counts at all! Could not fill energy matrices!");
    return;
  }

  std::vector<int> angle_bins(6,0);

  boost::multi_array<Real,6> reference(boost::extents[opts_.num_of_bins[0]][opts_.num_of_bins[1]][opts_.num_of_bins[2]][opts_.num_of_bins[3]][opts_.num_of_bins[4]][opts_.num_of_bins[5]]);

  if(reference_state=="classic"){
    for(int i=0;i<opts_.num_of_bins[0];++i){
      angle_bins[0]=i;
      for(int j=0;j<opts_.num_of_bins[1];++j){
        angle_bins[1]=j;
        for(int k=0;k<opts_.num_of_bins[2];++k){
          angle_bins[2]=k;
          for(int l=0;l<opts_.num_of_bins[3];++l){
            angle_bins[3]=l;
            for(int m=0;m<opts_.num_of_bins[4];++m){
              angle_bins[4]=m;
              for(int n=0;n<opts_.num_of_bins[5];++n){
                angle_bins[5]=n;
                reference[i][j][k][l][m][n] = stat->GetCount(angle_bins)/total_count;
              }
            }
          }
        }
      }
    }
  }
  else if(reference_state=="uniform"){
    Real n_o_b=1.0;
    for(int i=0;i<6;++i){
      n_o_b*=opts_.num_of_bins[i];
    }
    for(int i=0;i<opts_.num_of_bins[0];++i){
      angle_bins[0]=i;
      for(int j=0;j<opts_.num_of_bins[1];++j){
        angle_bins[1]=j;
        for(int k=0;k<opts_.num_of_bins[2];++k){
          angle_bins[2]=k;
          for(int l=0;l<opts_.num_of_bins[3];++l){
            angle_bins[3]=l;
            for(int m=0;m<opts_.num_of_bins[4];++m){
              angle_bins[4]=m;
              for(int n=0;n<opts_.num_of_bins[5];++n){
                angle_bins[5]=n;
                reference[i][j][k][l][m][n] = 1.0/n_o_b;
              }
            }
          }
        }
      }
    }
  }
  else{
    std::stringstream ss;
    ss << reference_state << " is not implemented as reference state in torsion potential!"
       << "implemented are: classic, uniform";
    throw io::IOException(ss.str());
  }

  Real e;
  Real propensity;

  for(std::vector<String>::iterator it=opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){
    Real total_count_group=stat->GetTotalCount(*it);
    for(int i=0;i<opts_.num_of_bins[0];++i){
      angle_bins[0]=i;
      for(int j=0;j<opts_.num_of_bins[1];++j){
        angle_bins[1]=j;
        for(int k=0;k<opts_.num_of_bins[2];++k){
          angle_bins[2]=k;
          for(int l=0;l<opts_.num_of_bins[3];++l){
            angle_bins[3]=l;
            for(int m=0;m<opts_.num_of_bins[4];++m){
              angle_bins[4]=m;
              for(int n=0;n<opts_.num_of_bins[5];++n){
                angle_bins[5]=n;
                propensity=0.0;
                if(reference[i][j][k][l][m][n]>0 && total_count_group>0){
                  propensity=stat->GetCount(*it,angle_bins)/total_count_group/reference[i][j][k][l][m][n];
                }
                e=log(1+opts_.sigma*total_count_group)-log(1+opts_.sigma*total_count_group*propensity);
                energies_[*it].Set(Index(i,j,k,l,m,n),e);
              }
            }
          }
        }
      }
    }
  }
}


TorsionEnergies* TorsionPotential::Data(std::vector<String>& residue_names){

  if(residue_names.size()!=3){
    std::stringstream ss;
    ss<<"requires 3 residue names, only ";
    ss<<residue_names.size();
    ss<<" provided.";
    throw std::runtime_error(ss.str());
  }

  String group_id=this->FindStat(residue_names);

  if(group_id == "") {
    String err = "Could not find Torsionenergies for provided residue names!";
    throw std::runtime_error(err);
  }

  std::map<String, TorsionEnergies>::iterator it = energies_.find(group_id);

  if(it == energies_.end()) {
    String err = "Could not find Torsionenergies for provided residue names!";
    throw std::runtime_error(err);
  }

  TorsionEnergies* ptr = &(it->second);

  return ptr;
}

}
