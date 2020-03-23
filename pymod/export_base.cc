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
#include <ost/mol/mol.hh>
#include <qmean/statistic_base.hh>
#include <qmean/potential_base.hh>
#include <vector>
#include <exception>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace ost;
using namespace ost::mol;
using namespace boost::python;

template<typename T1, typename T2>
T2 Get(T1 p, const String& key){
  std::vector<String> keys=p->GetKeys();
  if(std::find(keys.begin(), keys.end(), key)!=keys.end()){
    return (*p)[key];
  }
  std::stringstream ss;
  ss << "Key error! '" << key << "' does not exist in container!";
  throw std::runtime_error(ss.str());
}

template<typename T1, typename T2>
void Set(T1 p, const String key, T2 ptr){
  (*p)[key]=ptr;
}

template<typename T>
boost::python::list GetKeys(T p){
  std::vector<String> keys = p->GetKeys();
  boost::python::list ret_list = VecToList<String>(keys);
  return ret_list;
}

template<typename T>
void DeleteItem(T p, const String& key){
  p->Remove(key);
}

template<typename T>
bool ContainsItem(T p, const String& key){
  std::vector<String> keys = p->GetKeys();
  if(std::find(keys.begin(), keys.end(), key) != keys.end()) {
    return true;
  } 
  else {
    return false;
  }
} 

template<typename T>
int GetSize(T p){
  return p->GetSize();
}


void export_Base()
{

  class_<StatisticContainer>("StatisticContainer", init<>())
      .def("Load", &StatisticContainer::Load,(arg("filename"))).staticmethod("Load")
      .def("Save", &StatisticContainer::Save, (arg("filename")))
      .def("GetKeys", GetKeys<StatisticContainerPtr>)
      .def("__getitem__",Get<StatisticContainerPtr,StatisticBasePtr>)
      .def("__setitem__",Set<StatisticContainerPtr, StatisticBasePtr>)
      .def("__delitem__",DeleteItem<StatisticContainerPtr>)
      .def("__contains__",ContainsItem<StatisticContainerPtr>)
      .def("__len__",GetSize<StatisticContainerPtr>)

  ;

  class_<PotentialContainer>("PotentialContainer", init<>())
      .def("Load", &PotentialContainer::Load, (arg("filename"))).staticmethod("Load")
      .def("Save", &PotentialContainer::Save, (arg("filename")))
      .def("GetKeys", GetKeys<PotentialContainerPtr>)
      .def("__getitem__", Get<PotentialContainerPtr,PotentialBasePtr>)
      .def("__setitem__", Set<PotentialContainerPtr,PotentialBasePtr>)
      .def("__delitem__", DeleteItem<PotentialContainerPtr>)
      .def("__contains__",ContainsItem<PotentialContainerPtr>)
      .def("__len__",GetSize<PotentialContainerPtr>)
  ; 

  class_<StatisticBase, boost::noncopyable>("StatisticBase", no_init);

  class_<PotentialBase, boost::noncopyable>("PotentialBase", no_init);

  register_ptr_to_python<PotentialContainerPtr>();
  register_ptr_to_python<StatisticContainerPtr>();
  register_ptr_to_python<StatisticBasePtr>();
  register_ptr_to_python<PotentialBasePtr>();

}

