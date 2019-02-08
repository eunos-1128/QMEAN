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

