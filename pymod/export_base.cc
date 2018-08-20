// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

