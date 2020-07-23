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

#include <qmean/extract_data_helper.hh>

namespace qmean {

void ExtractTemplateDataDisCo(const String& entry_name, const String& chain_name,
                              const ost::seq::AlignmentHandle& aln,
                              ost::db::LinearIndexerPtr indexer,
                              ost::db::LinearCharacterContainerPtr seqres_container,
                              ost::db::LinearCharacterContainerPtr atomseq_container,
                              ost::db::LinearPositionContainerPtr& position_container,
                              std::vector<int>& residue_numbers,
                              geom::Vec3List& positions) {

  std::pair<uint64_t, uint64_t> data_range = indexer->GetDataRange(entry_name,
                                                                   chain_name);

  String template_seqres = aln.GetSequence(1).GetGaplessString();
  data_range.first += aln.GetSequence(1).GetOffset();
  data_range.second = data_range.first + template_seqres.size();

  // check, whether the the template seqres is consistent with what
  // we find in seqres_container
  String expected_template_seqres;
  seqres_container->GetCharacters(data_range, expected_template_seqres);
  if(expected_template_seqres != template_seqres) {
    throw std::runtime_error("Template sequence in input alignment is "
                             "inconsistent with sequence in SEQRES container!");
  }

  String template_atomseq;
  atomseq_container->GetCharacters(data_range, template_atomseq);
  geom::Vec3List extracted_positions;
  position_container->GetPositions(data_range, extracted_positions);

  // prepare output
  uint template_atomseq_size = template_atomseq.size();
  residue_numbers.clear();
  residue_numbers.reserve(template_atomseq_size);  
  positions.clear();
  positions.reserve(template_atomseq_size);

  // extract
  uint current_rnum = aln.GetSequence(0).GetOffset() + 1;
  uint current_template_pos = 0;
  String seqres_seq = aln.GetSequence(0).GetString();
  String template_seq = aln.GetSequence(1).GetString();

  for(int i = 0; i < aln.GetLength(); ++i) {
    if(seqres_seq[i] != '-' && template_seq[i] != '-') {
      if(template_atomseq[current_template_pos] != '-') {
        // it is aligned and we have a valid position!
        residue_numbers.push_back(current_rnum);
        positions.push_back(extracted_positions[current_template_pos]);
      }
    }

    if(seqres_seq[i] != '-') {
      ++current_rnum;
    }

    if(template_seq[i] != '-') {
      ++current_template_pos;
    }
  }
}

void ExtractTemplateDataGMQE(const String& entry_name, const String& chain_name,
                             const ost::seq::AlignmentHandle& aln,
                             ost::db::LinearIndexerPtr indexer, 
                             ost::db::LinearCharacterContainerPtr seqres_container,
                             ost::db::LinearCharacterContainerPtr atomseq_container,
                             ost::db::LinearCharacterContainerPtr dssp_container,
                             ost::db::LinearPositionContainerPtr& n_position_container,
                             ost::db::LinearPositionContainerPtr& ca_position_container,
                             ost::db::LinearPositionContainerPtr& c_position_container,
                             ost::db::LinearPositionContainerPtr& cb_position_container,
                             std::vector<int>& residue_numbers,
                             String& dssp,
                             geom::Vec3List& n_positions,
                             geom::Vec3List& ca_positions,
                             geom::Vec3List& c_positions,
                             geom::Vec3List& cb_positions) {

  std::pair<uint64_t, uint64_t> data_range = indexer->GetDataRange(entry_name,
                                                                   chain_name);

  String template_seqres = aln.GetSequence(1).GetGaplessString();
  data_range.first += aln.GetSequence(1).GetOffset();
  data_range.second = data_range.first + template_seqres.size();

  // check, whether the the template seqres is consistent with what
  // we find in seqres_container
  String expected_template_seqres;
  seqres_container->GetCharacters(data_range, expected_template_seqres);
  if(expected_template_seqres != template_seqres) {
    throw std::runtime_error("Template sequence in input alignment is "
                             "inconsistent with sequence in SEQRES container!");
  }

  String template_atomseq;
  atomseq_container->GetCharacters(data_range, template_atomseq);
  String extracted_dssp;
  dssp_container->GetCharacters(data_range, extracted_dssp);
  geom::Vec3List extracted_n_positions;
  n_position_container->GetPositions(data_range, extracted_n_positions);
  geom::Vec3List extracted_ca_positions;
  ca_position_container->GetPositions(data_range, extracted_ca_positions);
  geom::Vec3List extracted_c_positions;
  c_position_container->GetPositions(data_range, extracted_c_positions);
  geom::Vec3List extracted_cb_positions;
  cb_position_container->GetPositions(data_range, extracted_cb_positions);

  // prepare output
  uint template_atomseq_size = template_atomseq.size();
  residue_numbers.clear();
  residue_numbers.reserve(template_atomseq_size);  
  dssp.clear();
  dssp.reserve(template_atomseq_size);
  n_positions.clear();
  n_positions.reserve(template_atomseq_size);
  ca_positions.clear();
  ca_positions.reserve(template_atomseq_size);
  c_positions.clear();
  c_positions.reserve(template_atomseq_size);
  cb_positions.clear();
  cb_positions.reserve(template_atomseq_size);

  // extract
  uint current_rnum = aln.GetSequence(0).GetOffset() + 1;
  uint current_template_pos = 0;
  String seqres_seq = aln.GetSequence(0).GetString();
  String template_seq = aln.GetSequence(1).GetString();

  for(int i = 0; i < aln.GetLength(); ++i) {
    if(seqres_seq[i] != '-' && template_seq[i] != '-') {
      if(template_atomseq[current_template_pos] != '-') {
        // it is aligned and we have a valid position!
        residue_numbers.push_back(current_rnum);
        dssp.push_back(extracted_dssp[current_template_pos]);
        n_positions.push_back(extracted_n_positions[current_template_pos]);
        ca_positions.push_back(extracted_ca_positions[current_template_pos]);
        c_positions.push_back(extracted_c_positions[current_template_pos]);
        cb_positions.push_back(extracted_cb_positions[current_template_pos]);
      }
    }

    if(seqres_seq[i] != '-') {
      ++current_rnum;
    }

    if(template_seq[i] != '-') {
      ++current_template_pos;
    }
  }
}

} //ns

