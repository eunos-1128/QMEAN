// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
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
    .def_readonly("tilt", &FindMemParam::tilt)
    .def_readonly("angle", &FindMemParam::angle)
    .def_readonly("width", &FindMemParam::width)
    .def_readonly("pos", &FindMemParam::pos)
    .def_readonly("energy", &FindMemParam::energy)
    .def_readonly("axis", &FindMemParam::axis)
    .def_readonly("tilt_axis", &FindMemParam::tilt_axis)
    .def_readonly("membrane_axis", &FindMemParam::GetMembraneAxis)
    .def("Save", &FindMemParam::Save)
    .def("Load", &FindMemParam::Load).staticmethod("Load")
  ;

  class_<SolvationGrid>("SolvationGrid", init<geom::AlignedCuboid&, Real>())
    .def("GetOrigin", &SolvationGrid::GetOrigin)
    .def("AddView", &AddViewEntity, (arg("handle")))
    .def("AddView", &AddViewView, (arg("view")))
    .def("AddSurface", &SolvationGrid::AddSurface, (arg("surface")))
    .def("Flood", &SolvationGrid::Flood)
    .def("Flood3D", &SolvationGrid::Flood3D, (arg("start_pos")))
    .def("IsFilled", &SolvationGrid::IsFilled, (arg("pos")))
    //.def("AsIMG", &SolvationGrid::AsIMG)
    .def("GetSolvatedSurface", &SolvationGrid::GetSolvatedSurface)
    .def("GetSolvatedView", &SolvationGrid::GetSolvatedView)
  ;
  register_ptr_to_python<SolvationGridPtr>();

  def("FillMembraneDummies",&FillMembraneDummies,(arg("cuboid"),arg("surface"),arg("solvation_grid_bin_size")=0.5,arg("density")=0.035));
  def("FindMembrane",&FindMembraneWrapper,(arg("ent"),arg("surf"),arg("asa")));
  def("TransformSurface",&TransformSurface,(arg("surf"),arg("transform")));

}

