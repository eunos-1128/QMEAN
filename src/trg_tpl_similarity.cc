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

#include <ost/geom/vecmat3_op.hh>
#include <qmean/trg_tpl_similarity.hh>

namespace qmean {

TrgTplSimilarity::TrgTplSimilarity(const String& seqres,
                                   const std::vector<int>& residue_numbers, 
                                   const geom::Vec3List& ca_pos,
                                   Real distance_threshold) {

  if(ca_pos.size() != residue_numbers.size()) {
    throw std::runtime_error("ca_pos and residue_numbers inconsistent in size");
  }

  for(auto it = residue_numbers.begin(); it != residue_numbers.end(); ++it) {
    if(*it < 1 || *it > static_cast<int>(seqres.size())) {
      throw std::runtime_error("invalid residue number observed");
    }
  }

  ref_seq = seqres; 
  ref_distances.clear();
  ref_indices.clear();
  for(auto it = residue_numbers.begin(); it != residue_numbers.end(); ++it) {
    ref_indices.push_back((*it)-1);
  }
  n_interactions = 0;
  for(uint i = 0; i < ca_pos.size(); ++i) {
    std::vector<std::pair<int, Real> > distances;
    for(uint j = i + 1; j < ca_pos.size(); ++j) {
      Real d = geom::Distance(ca_pos[i], ca_pos[j]);
      if(d < distance_threshold) {
        distances.push_back(std::make_pair(ref_indices[j], d));
      }
    }
    ref_distances.push_back(distances);
    n_interactions += distances.size();
  }
}


Real TrgTplSimilarity::GetSimilarity(const std::vector<int>& residue_numbers,
                                     const geom::Vec3List& ca_pos) {

  if(ca_pos.size() != residue_numbers.size()) {
    throw std::runtime_error("ca_pos and residue_numbers inconsistent in size");
  }

  int max_resnum = ref_seq.size();
  for(auto it = residue_numbers.begin(); it != residue_numbers.end(); ++it) {
    if(*it < 1 || *it > max_resnum) {
      throw std::runtime_error("invalid residue number observed");
    }
  }

  std::vector<bool> valid_positions(ref_seq.size(), false);
  std::vector<geom::Vec3> positions(ref_seq.size(), geom::Vec3());
  for(uint i = 0; i < residue_numbers.size(); ++i) {
    valid_positions[residue_numbers[i]-1] = true;
    positions[residue_numbers[i]-1] = ca_pos[i];
  }

  int n_05 = 0;
  int n_1 = 0;
  int n_2 = 0;
  int n_4 = 0;
  for(uint i = 0; i < ref_distances.size(); ++i) {
    int ref_idx = ref_indices[i];
    if(valid_positions[ref_idx]) {
      for(auto it = ref_distances[i].begin(); it != ref_distances[i].end(); ++it) {
        if(valid_positions[it->first]) {
          Real d = geom::Distance(positions[ref_idx], positions[it->first]);
          Real d_diff = std::abs(d - it->second);
          n_4 += int((d_diff < Real(4.0)));
          n_2 += int((d_diff < Real(2.0)));
          n_1 += int((d_diff < Real(1.0)));
          n_05 += int((d_diff < Real(0.5)));
        }
      } 
    }
  }
  return Real(0.25) * (static_cast<Real>(n_4)/n_interactions + 
                       static_cast<Real>(n_2)/n_interactions +
                       static_cast<Real>(n_1)/n_interactions + 
                       static_cast<Real>(n_05)/n_interactions);
}

}

