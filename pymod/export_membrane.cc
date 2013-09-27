#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <qmean/vec_list_magic.hh>

#include <qmean/membrane.hh>


using namespace boost::python;


void FindMembraneWrapper(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, list& asa){
  std::vector<Real> real_asa = ListToVec<Real>(asa);
  FindMembrane(ent, surf, real_asa);
}






void export_Membrane()
{

  class_<SolvationGrid>("SolvationGrid", init<geom::AlignedCuboid&, Real>())
    .def("GetOrigin", &SolvationGrid::GetOrigin)
    .def("AddView", &SolvationGrid::AddView, (arg("view")))
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
}