// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
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

#include <qmean/gmqe_scores.hh>

namespace qmean {

GMQEScoreCalculator::GMQEScoreCalculator(PotentialContainerPtr potentials,
                                DisCoContainerPtr disco_container,
                                const ost::seq::SequenceHandle& seqres,
                                const ost::seq::SequenceHandle& psipred_pred,
                                const std::vector<int>& psipred_cfi) {

  // some consistency checks
  if(seqres.GetString() != disco_container->GetSeqres().GetString()) {
    throw std::runtime_error("Provided SEQRES does not match SEQRES in "
                             "DisCoContainer!");
  }

  // check whether the psipred_pred prediction only contains valid symbols 
  // [H, E, C] is further down...
  if(seqres.GetLength() != psipred_pred.GetLength()) {
    throw std::runtime_error("Provided SEQRES is inconsistent with provided " 
                             "psipred prediction!");
  }

  if(seqres.GetLength() != psipred_cfi.size()) {
    throw std::runtime_error("Provided SEQRES is inconsistent with provided " 
                             "psipred confidences!");    
  }

  // check, whether all required potentials are present
  std::vector<String> potential_keys = potentials->GetKeys();

  if(std::find(potential_keys.begin(), potential_keys.end(), "reduced_coil") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: reduced_coil to be "
                             "present in provided potential container!");
  }

  if(std::find(potential_keys.begin(), potential_keys.end(), "reduced_helix") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: reduced_helix to be "
                             "present in provided potential container!");
  }

  if(std::find(potential_keys.begin(), potential_keys.end(), "reduced_extended") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: reduced__extended to "
                             "be present in provided potential container!");
  }

  if(std::find(potential_keys.begin(), potential_keys.end(), "cb_packing_coil") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: cb_packing_coil to be "
                             "present in provided potential container!");
  }

  if(std::find(potential_keys.begin(), potential_keys.end(), "cb_packing_helix") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: cb_packing_helix to be "
                             "present in provided potential container!");
  }

  if(std::find(potential_keys.begin(), potential_keys.end(), "cb_packing_extended") ==
     potential_keys.end()) {
    throw std::runtime_error("Expect a potential with key: cb_packing_extended to "
                             "be present in provided potential container!");
  }
  
  // The potentials stored in the PotentialContainer are of basetype 
  // PotentialBase. Lets do some casting...
  qmean::ReducedPotentialPtr reduced_coil = 
  boost::dynamic_pointer_cast<qmean::ReducedPotential>(
  (*potentials)["reduced_coil"]);

  qmean::ReducedPotentialPtr reduced_helix = 
  boost::dynamic_pointer_cast<qmean::ReducedPotential>(
  (*potentials)["reduced_helix"]);

  qmean::ReducedPotentialPtr reduced_extended =
  boost::dynamic_pointer_cast<qmean::ReducedPotential>(
  (*potentials)["reduced_extended"]);

  qmean::CBPackingPotentialPtr cb_packing_coil = 
  boost::dynamic_pointer_cast<qmean::CBPackingPotential>(
  (*potentials)["cb_packing_coil"]);

  qmean::CBPackingPotentialPtr cb_packing_helix = 
  boost::dynamic_pointer_cast<qmean::CBPackingPotential>(
  (*potentials)["cb_packing_helix"]);

  qmean::CBPackingPotentialPtr cb_packing_extended =
  boost::dynamic_pointer_cast<qmean::CBPackingPotential>(
  (*potentials)["cb_packing_extended"]);  

  qmean::TorsionPotentialPtr torsion =
  boost::dynamic_pointer_cast<qmean::TorsionPotential>(
  (*potentials)["torsion"]);

  if(!(reduced_coil && reduced_helix && reduced_extended &&
       cb_packing_coil && cb_packing_helix && cb_packing_extended &&
       torsion)) {
    throw std::runtime_error("Could not apply dynamic cast to potential stored "
                             "in the provided PotentialContainer. They must be "
                             "of type CBPackingPotential,  "
                             "ReducedPotential and TorsionPotential!");
  }

  // check, whether the potential parametrizations are consistent
  qmean::impl::ReducedOpts reduced_coil_opts= reduced_coil->GetOpts();
  qmean::impl::ReducedOpts reduced_helix_opts = reduced_helix->GetOpts();
  qmean::impl::ReducedOpts reduced_extended_opts = reduced_extended->GetOpts();

  qmean::impl::CBPackingOpts cb_packing_coil_opts= cb_packing_coil->GetOpts();
  qmean::impl::CBPackingOpts cb_packing_helix_opts = cb_packing_helix->GetOpts();
  qmean::impl::CBPackingOpts cb_packing_extended_opts = cb_packing_extended->GetOpts();


  if(reduced_coil_opts != reduced_helix_opts || 
     reduced_coil_opts != reduced_extended_opts) {
    throw std::runtime_error("Options of provided reduced potentials must be "
                             " consistent!");
  }

  if(cb_packing_coil_opts != cb_packing_helix_opts || 
     cb_packing_coil_opts != cb_packing_extended_opts) {
    throw std::runtime_error("Options of provided cb_packing potentials must be "
                             " consistent!");
  }

  // to simplify the code later on, we assume a lower cutoff of 0.0 for the
  // reduced potential.
  if(reduced_coil_opts.lower_cutoff != 0.0) {
    throw std::runtime_error("The reduced potentials are assumed to have a "
                             "lower cutoff of 0.0!");
  }

  disco_container_ = disco_container;

  reduced_cutoff_ = reduced_coil_opts.upper_cutoff;
  reduced_squared_cutoff_ = reduced_cutoff_ * reduced_cutoff_;
  reduced_seq_sep_ = reduced_coil_opts.sequence_sep;
  cb_packing_squared_cutoff_ = cb_packing_coil_opts.cutoff * 
                               cb_packing_coil_opts.cutoff;
  cb_packing_max_count_ = cb_packing_coil_opts.max_counts;
  disco_cutoff_ = disco_container_->GetDistCutoff();
  disco_squared_cutoff_ = disco_cutoff_ * disco_cutoff_;
  disco_bin_size_ = disco_container_->GetBinSize();

  // We only store the raw data of the potentials. However, we keep a shared 
  // pointer as a member in order to keep the reference counts up and avoid 
  // invalidation of the raw data. 
  potential_container_ = potentials;

  // the actual amino acids defined in ost::conop is the only thing we need
  for(int i = 0; i < seqres.GetLength(); ++i) {
    ost::conop::AminoAcid aa = ost::conop::OneLetterCodeToAminoAcid(seqres[i]);
    if(aa == ost::conop::XXX) {
      throw std::runtime_error("Only standard amino acids allowed in SEQRES!"); 
    }
    amino_acids_.push_back(aa);
  }

  for(int i = 0; i < psipred_pred.GetLength(); ++i) {
    if(psipred_pred[i] == 'H') {
      cb_packing_energies_.push_back(cb_packing_helix->Data());
      reduced_energies_.push_back(reduced_helix->Data());
    }
    else if(psipred_pred[i] == 'E') {
      cb_packing_energies_.push_back(cb_packing_extended->Data());
      reduced_energies_.push_back(reduced_extended->Data());
    }
    else if(psipred_pred[i] == 'C') {
      cb_packing_energies_.push_back(cb_packing_coil->Data());
      reduced_energies_.push_back(reduced_coil->Data());
    }
    else {
      throw std::runtime_error("Expect only characters in ['H', 'E', 'C'] in "
                               "psipred_pred!");
    }
    psipred_pred_.push_back(psipred_pred[i]);
  }
  psipred_cfi_ = psipred_cfi;

  // the torsion potential is agnostic of secondary structure prediction...
  std::vector<String> aa_names;
  for(std::vector<ost::conop::AminoAcid>::iterator it = amino_acids_.begin();
      it != amino_acids_.end(); ++it) {
    aa_names.push_back(ost::conop::AminoAcidToResidueName(*it));
  }

  std::vector<String> aa_triplet;
  aa_triplet.push_back("ALA");
  aa_triplet.push_back("ALA");
  aa_triplet.push_back("ALA");

  // first will be invalid anyway...
  torsion_energies_.push_back(NULL);
  for(int i = 1; i < static_cast<int>(aa_names.size()) - 1; ++i) {
    aa_triplet[0] = aa_names[i-1];
    aa_triplet[1] = aa_names[i];
    aa_triplet[2] = aa_names[i+1];
    torsion_energies_.push_back(torsion->Data(aa_triplet));
  }
  // last one is again invalid...
  torsion_energies_.push_back(NULL);
}

void GMQEScoreCalculator::Eval(const geom::Vec3List& n_positions,
                               const geom::Vec3List& ca_positions,
                               const geom::Vec3List& c_positions,
                               const geom::Vec3List& cb_positions,
                               const ost::seq::SequenceHandle& dssp_states,
                               const std::vector<int>& residue_numbers,
                               Real& dist_const, Real& reduced, 
                               Real& cb_packing, Real& torsion, 
                               Real& ss_agreement, Real& coverage) const {

  // some consistency checks
  if(residue_numbers.size() != n_positions.size() ||
  	 residue_numbers.size() != ca_positions.size() ||
  	 residue_numbers.size() != c_positions.size() ||
  	 residue_numbers.size() != cb_positions.size()) {
    throw std::runtime_error("Expect resnums and positions to be of same size");
  }

  if(residue_numbers.size() != dssp_states.GetLength()) {
    throw std::runtime_error("Expect resnums and dssp states to be of same size!");
  }

  int seqres_size = amino_acids_.size();
  for(uint i = 0; i < residue_numbers.size(); ++i) {
  	if(residue_numbers[i] < 1 || residue_numbers[i] > seqres_size) {
  	  throw std::runtime_error("Invalid resnum");
  	}
  }

  int num_positions = ca_positions.size(); 
  int pos_idx = 0;
  std::vector<Real> squared_ca_distances(num_positions * num_positions / 2);
  for(int i = 0; i < num_positions; ++i) {
    for(int j = i + 1; j < num_positions; ++j, ++pos_idx) {
      squared_ca_distances[pos_idx] = geom::Length2(ca_positions[i] - 
                                                    ca_positions[j]);
    }
  }

  Real max_squared_ca_cutoff = std::max(reduced_squared_cutoff_,
                                        disco_squared_cutoff_);  

  std::vector<Real> reduced_scores(num_positions, 0.0);
  std::vector<Real> disco_scores(num_positions, 0.0);

  std::vector<int> cb_packing_counts(num_positions, 0);
  std::vector<int> reduced_counts(num_positions, 0);
  std::vector<int> disco_counts(num_positions, 0);

  geom::Vec3List direction_vectors(num_positions);
  for(int i = 0; i < num_positions; ++i) {
    direction_vectors[i] = 
    geom::Normalize(geom::Normalize(ca_positions[i]-n_positions[i]) +
                    geom::Normalize(ca_positions[i]-c_positions[i]));
  }

  pos_idx = 0;
  for(int i = 0; i < num_positions; ++i) {
    for(int j = i + 1; j < num_positions; ++j, ++pos_idx) {

      if(squared_ca_distances[pos_idx] < max_squared_ca_cutoff) {

        int rnum_i = residue_numbers[i];
        int rnum_j = residue_numbers[j];

        // let the scoring action begin
        Real ca_distance = std::sqrt(squared_ca_distances[pos_idx]);
        Real squared_cb_distance = geom::Length2(cb_positions[i] - 
                                                 cb_positions[j]);

        if(squared_cb_distance < cb_packing_squared_cutoff_) {
          cb_packing_counts[i] += 1;
          cb_packing_counts[j] += 1;
        }

        if(ca_distance < disco_cutoff_ && 
           disco_container_->HasConstraint(rnum_i, rnum_j)) {

          const std::vector<Real>& constraint = 
          disco_container_->GetConstraint(rnum_i, rnum_j);

          // we perform a linear interpolation
          uint bin_lower = static_cast<uint>(ca_distance / disco_bin_size_);
          uint bin_upper = bin_lower + 1; 
          Real w_one = (bin_upper * disco_bin_size_ - ca_distance) / 
                       disco_bin_size_;
          Real w_two = Real(1.0) - w_one;
          Real score = w_one * constraint[bin_lower] + 
                       w_two * constraint[bin_upper];

          disco_scores[i] += score;
          disco_scores[j] += score;
          disco_counts[i] += 1;
          disco_counts[j] += 1;
        }

        if(ca_distance < reduced_cutoff_ && 
           std::abs(rnum_i - rnum_j) >= reduced_seq_sep_) {

          geom::Vec3 v_i = direction_vectors[i];
          geom::Vec3 v_j = direction_vectors[j];
          geom::Vec3 v_ij = (ca_positions[j] - ca_positions[i]);
          v_ij /= ca_distance;

          Real dot_product = geom::Dot(v_i, v_ij);
          dot_product = std::max(Real(-0.999999), dot_product);
          dot_product = std::min(Real(1.0), dot_product);
          Real alpha = std::acos(dot_product);

          dot_product = -geom::Dot(v_ij, v_j);
          dot_product = std::max(Real(-0.999999), dot_product);
          dot_product = std::min(Real(1.0), dot_product);
          Real beta = std::acos(dot_product);

          geom::Vec3 v_i_ij_cross = -geom::Cross(v_i, v_ij);
          geom::Vec3 v_ij_j_cross = geom::Cross(v_ij, v_j);
          Real gamma = std::atan2(-geom::Dot(v_i, v_ij_j_cross),
                                  geom::Dot(v_i_ij_cross,v_ij_j_cross));

          ost::conop::AminoAcid aa_i = amino_acids_[rnum_i - 1];
          ost::conop::AminoAcid aa_j = amino_acids_[rnum_j - 1];
          
          Real e_i = reduced_energies_[rnum_i-1]->Get(aa_i, aa_j, ca_distance, 
                                                      alpha, beta, gamma);
          Real e_j = reduced_energies_[rnum_j-1]->Get(aa_j, aa_i, ca_distance, 
                                                      beta, alpha, gamma);

          reduced_scores[i] += e_i;
          reduced_scores[j] += e_j;

          reduced_counts[i] += 1;
          reduced_counts[j] += 1;
        }
      }
    }
  }

  Real summed_reduced = 0.0;
  Real summed_cb_packing = 0.0;
  Real summed_disco = 0.0;
  Real summed_ss_agreement = 0.0;

  for(int i = 0; i < num_positions; ++i) {

    // for cb_packing we still need to extract the actual energies
    int rnum = residue_numbers[i];
    ost::conop::AminoAcid aa = amino_acids_[rnum - 1];
    int c = std::min(cb_packing_counts[i], cb_packing_max_count_); 
    summed_cb_packing += cb_packing_energies_[rnum-1]->Get(aa, c);

    summed_ss_agreement += ss_agreement_.LogOddScore(dssp_states[i],
                                                     psipred_pred_[rnum-1],
                                                     psipred_cfi_[rnum-1]);
    
    // the other two are trivial
    if(reduced_counts[i] > 0) {
      summed_reduced += reduced_scores[i] / reduced_counts[i];  
    }

    if(disco_counts[i] > 0) {
      summed_disco += disco_scores[i] / disco_counts[i];
    }
  }  

  // torsion is done separately
  std::vector<Real> 
  phi_angles(num_positions, std::numeric_limits<Real>::quiet_NaN());
  std::vector<Real> 
  psi_angles(num_positions, std::numeric_limits<Real>::quiet_NaN());

  // do phi angles
  for(int i = 1; i < num_positions; ++i) {
    if(residue_numbers[i-1]+1 == residue_numbers[i] &&
       geom::Length2(c_positions[i-1]-n_positions[i]) < Real(9.0)) {
      phi_angles[i] = geom::DihedralAngle(c_positions[i-1], n_positions[i],
                                          ca_positions[i], c_positions[i]);
    }
  }

  // do psi angles
  for(int i = 0; i < num_positions-1; ++i) {
    if(residue_numbers[i]+1 == residue_numbers[i+1] &&
       geom::Length2(c_positions[i]-n_positions[i+1]) < Real(9.0)) {
      psi_angles[i] = geom::DihedralAngle(n_positions[i], ca_positions[i],
                                          c_positions[i], n_positions[i+1]);
    }
  }

  Real summed_torsion = 0.0;
  int n_torsions = 0;
  for(int i = 1; i < num_positions-1; ++i) {
    if(phi_angles[i-1] == phi_angles[i-1] &&
       psi_angles[i-1] == psi_angles[i-1] &&
       phi_angles[i] == phi_angles[i] &&
       psi_angles[i] == psi_angles[i] &&
       phi_angles[i+1] == phi_angles[i+1] &&
       psi_angles[i+1] == psi_angles[i+1]) {
      summed_torsion += 
      torsion_energies_[residue_numbers[i]-1]->Get(phi_angles[i-1], 
                                                   psi_angles[i-1],
                                                   phi_angles[i], 
                                                   psi_angles[i],
                                                   phi_angles[i+1], 
                                                   psi_angles[i+1]);
      ++n_torsions;
    } 
  }
  
  dist_const = summed_disco;
  reduced = summed_reduced;
  cb_packing = summed_cb_packing;
  torsion = summed_torsion;
  ss_agreement = summed_ss_agreement;

  if(num_positions > 0) {
    dist_const /= num_positions;
    reduced /= num_positions;
    cb_packing /= num_positions;
    ss_agreement /= num_positions;
  }

  if(n_torsions > 0) {
    torsion /= n_torsions;
  }

  coverage = static_cast<Real>(num_positions) / amino_acids_.size();
}


} // ns
