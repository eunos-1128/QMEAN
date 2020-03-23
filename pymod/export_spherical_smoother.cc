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

