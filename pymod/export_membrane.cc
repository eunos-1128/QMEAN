#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <qmean/vec_list_magic.hh>

#include <qmean/membrane.hh>


using namespace boost::python;


FindMemParam FindMembraneWrapper(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, list& asa){
  std::vector<Real> real_asa = ListToVec<Real>(asa);
  return FindMembrane(ent, surf, real_asa);
}

void AddViewEntity(SolvationGrid& g, ost::mol::EntityHandle& e) { g.AddView(e); }
void AddViewView(SolvationGrid& g, ost::mol::EntityView& v) { g.AddView(v); }





void export_Membrane()
{

  class_<FindMemParam>("FindMemParam", no_init)
    .def_readwrite("tilt", &FindMemParam::tilt)
    .def_readwrite("angle", &FindMemParam::angle)
    .def_readwrite("width", &FindMemParam::width)
    .def_readwrite("pos", &FindMemParam::pos)
    .def_readwrite("euler_one", &FindMemParam::euler_one)
    .def_readwrite("euler_two", &FindMemParam::euler_two)
    .def_readwrite("energy", &FindMemParam::energy)
    .def_readonly("membrane_axis", &FindMemParam::GetMembraneAxis)
    .def_readonly("tilt_axis", &FindMemParam::GetTiltAxis)
  ;

  class_<SolvationGrid>("SolvationGrid", init<geom::AlignedCuboid&, Real>())
    .def("GetOrigin", &SolvationGrid::GetOrigin)
    .def("AddView", &AddViewEntity, (arg("handle")))
    .def("AddView", &AddViewView, (arg("view")))
    .def("AddSurface", &SolvationGrid::AddSurface, (arg("surface")))
    .def("Flood", &SolvationGrid::Flood)
    .def("IsFilled", &SolvationGrid::IsFilled, (arg("pos")))
    .def("AsIMG", &SolvationGrid::AsIMG)
    .def("GetSolvatedSurface", &SolvationGrid::GetSolvatedSurface)
    .def("GetSolvatedView", &SolvationGrid::GetSolvatedView)
  ;
  register_ptr_to_python<SolvationGridPtr>();

  def("FillMembraneDummies",&FillMembraneDummies,(arg("cuboid"),arg("surface"),arg("solvation_grid_bin_size")=0.5,arg("density")=0.035));
  def("FindMembrane",&FindMembraneWrapper,(arg("ent"),arg("surf"),arg("asa")));
  def("TransformSurface",&TransformSurface,(arg("surf"),arg("transform")));

}