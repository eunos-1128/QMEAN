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

