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

#ifndef QMEAN_EXTRACT_DATA_HELPER_HH
#define QMEAN_EXTRACT_DATA_HELPER_HH

#include <ost/geom/vec3.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/db/linear_indexer.hh>
#include <ost/db/binary_container.hh>
#include <ost/seq/alignment_handle.hh>

namespace qmean {

void ExtractTemplateDataDisCo(const String& entry_name, const String& chain_name,
                              const ost::seq::AlignmentHandle& aln,
                              ost::db::LinearIndexerPtr indexer, 
                              ost::db::LinearCharacterContainerPtr seqres_container,
                              ost::db::LinearCharacterContainerPtr atomseq_container,
                              ost::db::LinearPositionContainerPtr& position_container,
                              std::vector<int>& residue_numbers,
                              geom::Vec3List& positions);

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
                             geom::Vec3List& cb_positions);
} //ns

#endif
