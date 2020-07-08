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

#ifndef TRG_TPL_SIMILARITY_HH
#define TRG_TPL_SIMILARITY_HH

#include <vector>
#include <ost/geom/vec3.hh>



namespace qmean {

struct TrgTplSimilarity;
typedef boost::shared_ptr<TrgTplSimilarity> TrgTplSimilarityPtr;


struct TrgTplSimilarity{

TrgTplSimilarity(const String& seqres,
                 const std::vector<int>& residue_numbers, 
                 const geom::Vec3List& ca_pos,
                 Real distance_threshold = 15.0);

Real GetSimilarity(const std::vector<int>& residue_numbers,
                   const geom::Vec3List& ca_pos);

String ref_seq;
std::vector<std::vector<std::pair<int, Real> > > ref_distances;
std::vector<int> ref_indices;
int n_interactions;
};


}

#endif

