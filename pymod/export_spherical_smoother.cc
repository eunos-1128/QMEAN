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
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/spherical_smoother.hh>
#include <qmean/vec_list_magic.hh>

using namespace boost::python;
using namespace qmean;

boost::python::list WrapSmooth(SphericalSmootherPtr s, boost::python::list& v){
  std::vector<Real> values = ListToVec<Real>(v);
  std::vector<Real> result = s->Smooth(values);
  return VecToList<Real>(result);
}

SphericalSmootherPtr WrapSphericalSmootherConstructor(boost::python::list p, Real std){
  std::vector<geom::Vec3> pos = ListToVec<geom::Vec3>(p);
  SphericalSmootherPtr ss(new SphericalSmoother(pos, std));
  return ss;
}


void export_SphericalSmoother()
{ 
  class_<SphericalSmoother>("SphericalSmoother")
    .def("__init__",(arg("positions"),arg("std"),make_constructor(&WrapSphericalSmootherConstructor)))
    .def("Smooth", &WrapSmooth, (arg("values")))
  ;
  register_ptr_to_python<SphericalSmootherPtr>();
}

