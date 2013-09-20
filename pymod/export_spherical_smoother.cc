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