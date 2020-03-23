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

#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/clash_score.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;



namespace {
list WrapGetClashScores(ost::mol::EntityView& target, ost::mol::EntityView& env){
  std::vector<Real> scores = GetClashScores(target, env);
  list ret_list=VecToList<Real>(scores);
  return ret_list;
}
}

void export_clash() {
  def("GetClashScores", &GetClashScores, (arg("target"), arg("env")));
}
