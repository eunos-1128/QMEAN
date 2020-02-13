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

#ifndef QMEAN_INDEX_HH
#define QMEAN_INDEX_HH

/*
  Author: Marco Biasini
 */ 

#include <ost/stdint.hh>

#include <cstring>
#include <boost/shared_ptr.hpp>
#include <qmean/module_config.hh>

namespace qmean {

namespace impl {
template <uint32_t D>
class IndexBase {
public:
  enum { Dimension = D };
   IndexBase(const IndexBase& rhs) {
    memcpy(data_, rhs.data_, sizeof(uint32_t[D]));
  }
  IndexBase() {
    memset(data_, 0, sizeof(uint32_t[D]));
  }
  IndexBase& operator=(const IndexBase& rhs) {
    memcpy(data_, rhs.data_, sizeof(uint32_t[D]));
    return *this;
  }
  uint32_t operator[](uint32_t idx) const {
    assert(idx <= D);
    return data_[idx];
  }
  uint32_t& operator[](uint32_t idx) {
    assert(idx <= D);
    return data_[idx];
  }  
private:
  uint32_t data_[D];   
};

} // namespace impl

template <uint32_t D>
class Index;

template <>
class Index<1> : public impl::IndexBase<1> {
public:
  Index() : impl::IndexBase<1>() {}
  Index(uint32_t a) {
    (*this)[0]=a;
  }
};
template <>
class Index<2> : public impl::IndexBase<2> {
public:
  Index() : impl::IndexBase<2>() {}
  Index(uint32_t a, uint32_t b) {
    (*this)[0]=a;
    (*this)[1]=b;    
  }
};

template <>
class Index<3> : public impl::IndexBase<3> {
public:
  Index() : impl::IndexBase<3>() {}
  Index(uint32_t a, uint32_t b, uint32_t c) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
  }
};

template <>
class Index<4> : public impl::IndexBase<4> {
public:
  Index() : impl::IndexBase<4>() {}
  Index(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
  }
};

template <>
class Index<5> : public impl::IndexBase<5> {
public:
  Index() : impl::IndexBase<5>() {}
  Index(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;                
  }
};

template <>
class Index<6> : public impl::IndexBase<6> {
public:
  Index() : impl::IndexBase<6>() {}
  Index(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e, uint32_t f) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;
    (*this)[5]=f;    
  }
};

template <>
class Index<7> : public impl::IndexBase<7> {
public:
  Index() : impl::IndexBase<7>() {}
  Index(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e, uint32_t f, uint32_t g) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;
    (*this)[5]=f;    
    (*this)[6]=g;        
  }
};

template<uint32_t D>
class IndexIterator {
public:
  typedef Index<D>  IndexType;
  IndexIterator(const IndexType& s, const IndexType& e)
    : start_(s), end_(e), current_(s) {
    
  }  
  
  IndexIterator<D>& operator++() {
    uint32_t current_it=0;
    while (++current_[current_it] > end_[current_it]) {
      current_it++;
      if (current_it < D) {
        current_[current_it-1] = start_[current_it-1];
      } else {
        break;
      }
    }
    return *this;
  }
  const IndexType& operator *() const {
    return current_;
  }
  bool AtEnd() {
    return current_[D-1] > end_[D-1];
  }
private:
  IndexType  start_;
  IndexType  end_;
  IndexType  current_;

};

}

#endif

