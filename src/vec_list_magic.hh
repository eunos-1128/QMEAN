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

#ifndef VEC_LIST_MAGIC_HH
#define VEC_LIST_MAGIC_HH

#include <boost/python.hpp>
#include <vector>


template<typename T>
boost::python::list VecToList(const std::vector<T>& vec){
  boost::python::list l;
  for(typename std::vector<T>::const_iterator it=vec.begin();
      it!=vec.end(); ++it){
    l.append(*it);
  }
  return l;
}

template<typename T>
std::vector<T> ListToVec(const boost::python::list& l){

  std::vector<T> vec;
  for (int i = 0; i < boost::python::len(l); ++i){
    vec.push_back(boost::python::extract<T>(l[i]));
  }
  return vec;
}

#endif

