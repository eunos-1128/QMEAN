#include <qmean/hbond_impl.hh>

namespace qmean{ namespace impl{

void HBondPotentialImpl::SetEnvironment(const ost::mol::EntityView& env){

  env_.Clear();
  ost::mol::ResidueViewList res_list = env.GetResidueList();

  std::vector<int> states(res_list.size());
  for(uint i = 0; i < res_list.size(); ++i){
    states[i] = this->GetState(res_list[i].GetHandle());
  }

  //per definition, for being part of a particular state, at least one
  //neighbouring residue has to be of the same state
  
  //first do the first and last residue
  if(res_list.size() > 2){
    if(states[0] == 1){
      if(states[1] != 1) states[0] = 0;
    }
    if(states[0] == 2){
      if(states[1] != 2) states[0] = 0;    
    }
    if(states.back() == 1){
      if(states[states.size()-2] != 1) states.back() = 0;
    }
    if(states.back() == 2){
      if(states[states.size()-2] != 2) states.back() = 0;
    }
  }

  //let's do the remaining residues in between
  for(uint i = 1; i < states.size()-1; ++i){
    if(states[i] == 0) continue;
    if(states[i] == 1){
      if(states[i-1] == 1 || states[i+1] == 1) continue;
    }
    if(states[i] == 2){
      if(states[i-1] == 2 || states[i+1] == 2) continue;
    }
    states[i] = 0;
  }

  ost::mol::AtomView n,ca,c,o;
  geom::Vec3 n_pos,ca_pos,c_pos,o_pos;
  geom::Vec3 h_pos;
  Real phi;
  ost::mol::TorsionHandle phi_torsion;
  Real h_bond_length = 0.9970;
  Real h_angle = 2.0420;
  String chain_name;
  int num;

  for(uint i = 0; i < res_list.size(); ++i){
    n = res_list[i].FindAtom("N");
    ca = res_list[i].FindAtom("CA");
    c = res_list[i].FindAtom("C");
    o = res_list[i].FindAtom("O");

    if(!(n.IsValid() && ca.IsValid() && c.IsValid() && o.IsValid())) continue;

    n_pos = n.GetPos();
    ca_pos = ca.GetPos();
    c_pos = c.GetPos();
    o_pos = o.GetPos();
    phi = -1.0472;
    phi_torsion = res_list[i].GetHandle().GetPhiTorsion();
    if(phi_torsion.IsValid()) phi = phi_torsion.GetAngle();
    chain_name = res_list[i].GetChain().GetName();
    num = res_list[i].GetNumber().GetNum();
    this->ConstructHPos(h_bond_length,h_angle,phi+M_PI,
                        c_pos, ca_pos, n_pos, h_pos);

    HBondSpatialOrganizerItem new_item(n_pos,ca_pos, c_pos,
                                       o_pos,h_pos,states[i],chain_name,num);
    env_.Add(new_item,ca_pos);
  }
}

bool HBondPotentialImpl::VisitResidue(const ost::mol::ResidueHandle& res){

  ost::mol::AtomHandle n = res.FindAtom("N");
  ost::mol::AtomHandle ca = res.FindAtom("CA");
  ost::mol::AtomHandle c = res.FindAtom("C");
  ost::mol::AtomHandle o = res.FindAtom("O");

  if(!(n.IsValid() && ca.IsValid() && c.IsValid() && o.IsValid())) return false;

  geom::Vec3 h_pos;
  geom::Vec3 n_pos = n.GetPos();
  geom::Vec3 ca_pos = ca.GetPos();
  geom::Vec3 c_pos = c.GetPos();
  geom::Vec3 o_pos = o.GetPos();

  Real phi = -1.0472;
  ost::mol::TorsionHandle phi_torsion = res.GetPhiTorsion();
  if(phi_torsion.IsValid()) phi = phi_torsion.GetAngle();
  Real h_bond_length = 0.9970;
  Real h_angle = 2.0420;
  this->ConstructHPos(h_bond_length,h_angle,phi+M_PI,
                      c_pos, ca_pos, n_pos, h_pos);


  int state = this->GetState(res);
  ost::mol::ResidueHandle prev = res.GetPrev();
  ost::mol::ResidueHandle next = res.GetNext();
  int prev_state = 0;
  int next_state = 0;

  if(prev.IsValid()){
    prev_state = this->GetState(prev);
  }
  if(next.IsValid()){
    next_state = this->GetState(next);
  }

  switch(state){
    case 1:{
      if(prev_state != 1 && next_state != 1) state = 0;
      break;
    }
    case 2:{
      if(prev_state != 2 && next_state != 2) state = 0;
      break;      
    }
    default: break;
  }

  Real d,alpha,beta,gamma;
  int pairwise_state;
  std::vector<HBondSpatialOrganizerItem> in_reach = env_.FindWithin(ca_pos, 8.0);
  int num = res.GetNumber().GetNum();
  String chain_name = res.GetChain().GetName();
  int diff;

  //current residue as donor
  for(std::vector<HBondSpatialOrganizerItem>::iterator i = in_reach.begin();
      i != in_reach.end(); ++i){

    if(unique_bonds_){
      if(i->num-num<opts_.seq_sep){
        if(chain_name==i->chain_name){
          continue;
        }
      }
    }
    else{
      if(std::abs(i->num-num)<opts_.seq_sep){
        if(chain_name==i->chain_name){
          continue;
        }
      }
    }

    d = geom::Distance(h_pos,i->o_pos);
    if(d<opts_.d_min || d>=opts_.d_max) continue;
    alpha = geom::Angle(n_pos-h_pos,i->o_pos-h_pos);
    if(alpha<opts_.alpha_min || alpha>=opts_.alpha_max) continue;
    beta = geom::Angle(h_pos-i->o_pos,i->c_pos-i->o_pos);
    if(beta<opts_.beta_min || beta>=opts_.beta_max) continue;
    gamma = geom::DihedralAngle(i->ca_pos,i->c_pos,i->o_pos,h_pos);
    if(gamma<opts_.gamma_min || gamma>=opts_.gamma_max) continue;
    pairwise_state = 0;
    if(state == i->state) pairwise_state = state;

    this->OnInteraction(pairwise_state,d,alpha,beta,gamma);
  }

  //current residue as acceptor
  for(std::vector<HBondSpatialOrganizerItem>::iterator i = in_reach.begin();
      i != in_reach.end(); ++i){

    if(unique_bonds_){
      if(i->num-num<opts_.seq_sep){
        if(chain_name==i->chain_name){
          continue;
        }
      }
    }
    else{
      if(std::abs(i->num-num)<opts_.seq_sep){
        if(chain_name==i->chain_name){
          continue;
        }
      }
    }

    d = geom::Distance(o_pos,i->h_pos);
    if(d<opts_.d_min || d>=opts_.d_max) continue;
    alpha = geom::Angle(o_pos-i->h_pos,i->n_pos-i->h_pos);
    if(alpha<opts_.alpha_min || alpha>=opts_.alpha_max) continue;
    beta = geom::Angle(c_pos-o_pos,i->h_pos-o_pos);
    if(beta<opts_.beta_min || beta>=opts_.beta_max) continue;
    gamma = geom::DihedralAngle(ca_pos,c_pos,o_pos,i->h_pos);
    if(gamma<opts_.gamma_min || gamma>=opts_.gamma_max) continue;
    pairwise_state = 0;
    if(state == i->state) pairwise_state = state;
    this->OnInteraction(pairwise_state,d,alpha,beta,gamma);
  }

  return false;
}

int HBondPotentialImpl::GetState(const ost::mol::ResidueHandle& res){

  Real phi = -1.0472;
  Real psi = -0.78540;

  ost::mol::TorsionHandle phi_torsion = res.GetPhiTorsion();
  ost::mol::TorsionHandle psi_torsion = res.GetPsiTorsion();

  if(phi_torsion.IsValid()) phi = phi_torsion.GetAngle();
  if(psi_torsion.IsValid()) psi = psi_torsion.GetAngle();

  //check whether phi/psi pair is typical for
  //helix => -180 < phi < -20 ; -90 < psi < -10
  if(-M_PI < phi && phi < -0.34907 &&
     -M_PI/2 < psi && psi < -0.17453) return 1;

  //check whether phi/psi pair is typical for
  //sheet => -180 < phi < -20 ; 20 < psi < 180 or -180 < psi < -170
  if(-M_PI < phi && phi < -0.34907 &&
     (( 0.34907 < psi && psi < M_PI) || 
      (-M_PI < psi && psi < -2.9671))) return 2; 

  return 0;
}

void HBondPotentialImpl::ConstructHPos(Real bond_length, Real angle, Real dihedral, 
                                       geom::Vec3& A, geom::Vec3& B, geom::Vec3& C,
                                       geom::Vec3& pos){

  Real x,xx,x1;
  x = std::tan((M_PI-angle)/2);
  xx = x*x;
  x1 = 1+xx;
  Real bond_length_x_sin_angle = bond_length*2*x/x1;
  Real bond_length_x_cos_angle = bond_length*(1-xx)/x1;

  x = std::tan(dihedral/2);
  xx = x*x;
  x1 = 1+xx;
  Real sin_dihedral = 2*x/x1;
  Real cos_dihedral = (1-xx)/x1;

  Real ab[] = {B[0]-A[0],B[1]-A[1],B[2]-A[2]};

  Real norm_bc[] = {C[0]-B[0],C[1]-B[1],C[2]-B[2]};
  Real a = norm_bc[0] * norm_bc[0];
  a += norm_bc[1] * norm_bc[1];
  a += norm_bc[2] * norm_bc[2];
  a = 1.0 / std::sqrt(a);

  norm_bc[0]*=a;
  norm_bc[1]*=a;
  norm_bc[2]*=a;

  Real norm_n[] = {ab[1]*norm_bc[2]-norm_bc[1]*ab[2],
                   ab[2]*norm_bc[0]-norm_bc[2]*ab[0],
                   ab[0]*norm_bc[1]-norm_bc[0]*ab[1]};

  a = norm_n[0]*norm_n[0];  
  a += norm_n[1]*norm_n[1];  
  a += norm_n[2]*norm_n[2];  
  a = 1.0 / std::sqrt(a);

  norm_n[0]*=a;
  norm_n[1]*=a;
  norm_n[2]*=a;

  Real n_x_bc[] = {norm_n[1]*norm_bc[2]-norm_bc[1]*norm_n[2],
                   norm_n[2]*norm_bc[0]-norm_bc[2]*norm_n[0],
                   norm_n[0]*norm_bc[1]-norm_bc[0]*norm_n[1]};

  pos[0] = norm_bc[0]*bond_length_x_cos_angle + 
           n_x_bc[0]*cos_dihedral*bond_length_x_sin_angle + 
           norm_n[0]*sin_dihedral*bond_length_x_sin_angle + C[0];

  pos[1] = norm_bc[1]*bond_length_x_cos_angle + 
           n_x_bc[1]*cos_dihedral*bond_length_x_sin_angle + 
           norm_n[1]*sin_dihedral*bond_length_x_sin_angle + C[1];

  pos[2] = norm_bc[2]*bond_length_x_cos_angle + 
           n_x_bc[2]*cos_dihedral*bond_length_x_sin_angle + 
           norm_n[2]*sin_dihedral*bond_length_x_sin_angle + C[2];

}

}}
