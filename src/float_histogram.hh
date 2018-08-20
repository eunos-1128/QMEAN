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

#ifndef QMEAN_FLOAT_HISTOGRAM_HH
#define QMEAN_FLOAT_HISTOGRAM_HH
/*
  Author: Marco Biasini
 */
#include <ost/stdint.hh>

#include <qmean/module_config.hh> 
#include <qmean/multi_classifier.hh>

namespace qmean {

/// \brief histogram specialization for multi classifier  
template <typename T1,
          typename T2=impl::NullType,
          typename T3=impl::NullType,
          typename T4=impl::NullType,
          typename T5=impl::NullType,
          typename T6=impl::NullType,
          typename T7=impl::NullType>
class FloatHistogram : public MultiClassifier<Real, T1, T2, T3, T4, T5, T6, T7> 
{
public:
  typedef Classifier<T1>                       C1;
  typedef Classifier<T2>                       C2;
  typedef Classifier<T3>                       C3;
  typedef Classifier<T4>                       C4;
  typedef Classifier<T5>                       C5;
  typedef Classifier<T6>                       C6;
  typedef Classifier<T7>                       C7;  
#ifdef _MSC_VER
  FloatHistogram(typename C1::ConstRefType c1,
            typename C2::ConstRefType c2=C2::Type(),
            typename C3::ConstRefType c3=C3::Type(),
            typename C4::ConstRefType c4=C4::Type(),
            typename C5::ConstRefType c5=C5::Type(),
            typename C6::ConstRefType c6=C6::Type(),
            typename C7::ConstRefType c7=C7::Type())
#else
  FloatHistogram(typename C1::ConstRefType c1,
            typename C2::ConstRefType c2=typename C2::Type(),
            typename C3::ConstRefType c3=typename C3::Type(),
            typename C4::ConstRefType c4=typename C4::Type(),
            typename C5::ConstRefType c5=typename C5::Type(),
            typename C6::ConstRefType c6=typename C6::Type(),
            typename C7::ConstRefType c7=typename C7::Type())
#endif
   : MultiClassifier<Real, T1, T2, T3, T4, 
                     T5, T6, T7>(0, c1, c2, c3, c4, c5, c6, c7) {        
  }
  
  FloatHistogram(): 
    MultiClassifier<Real, T1, T2, T3, T4, T5, T6, T7>() 
  {
  }
};

} // ost::qa

#endif

