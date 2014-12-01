#include <ost/mol/mol.hh>

#include <ost/geom/geom.hh>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>
#include <limits>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <qmean/module_config.hh>
#include <ost/mol/surface_handle.hh>
//#include <ost/img/image_handle.hh>
#include <math.h>
#include <ost/img/image_factory.hh>
#include <ost/img/point.hh>
#include <ost/img/size.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/geom/transform.hh>
#include <ost/io/binary_data_source.hh>
#include <ost/io/binary_data_sink.hh>
#include <ost/io/io_exception.hh>

#include <ost/img/alg/levenberg_marquardt.h>


#include <ost/io/pdb_writer.hh>
#include <ost/io/io_profile.hh>


#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/SVD>
#include <Eigen/LU>

#include <boost/filesystem.hpp>


class SolvationGrid;
typedef boost::shared_ptr<SolvationGrid> SolvationGridPtr;

struct FindMemParam{
  FindMemParam() { }

  geom::Vec3 GetMembraneAxis();
  geom::Vec3 axis;
  geom::Vec3 tilt_axis;
  Real tilt;
  Real angle;
  Real width;
  Real pos;
  Real energy;

  static FindMemParam Load(const String& filename);
  void Save(const String& filename);

  template <typename DS>
  void Serialize(DS& ds){
    ds & tilt;
    ds & angle;
    ds & width;
    ds & pos;
    ds & energy;
    ds & axis[0];
    ds & axis[1];
    ds & axis[2];
    ds & tilt_axis[0];
    ds & tilt_axis[1];
    ds & tilt_axis[2];
  }




};

class DLLEXPORT_QMEAN SolvationGrid {
public:
  SolvationGrid() {};

  SolvationGrid(geom::AlignedCuboid cub, Real bin_size);

  void AddView(const ost::mol::EntityView& view) { view_ = view; };

  void AddView(const ost::mol::EntityHandle& handle) { this->AddView(handle.CreateFullView()); }

  void AddSurface(const ost::mol::SurfaceHandle& surface);

  void Flood();

  bool IsFilled(geom::Vec3& pos);

  geom::Vec3 GetOrigin() { return origin_; }

  //ost::img::ImageHandle AsIMG();

  void CalculateSolvatedSurface();

  void CalculateSolvatedView();

  ost::mol::SurfaceHandle GetSolvatedSurface();

  ost::mol::EntityView GetSolvatedView();

  ost::mol::SurfaceVertexIDList GetSolvatedVertices() { return solvated_vertices_; }


private:

  void FloodLevel(int level, std::pair<int, int> start_coordinates, int orig_value, int dest_value);

  void FloodLevel8Neighbours(int level, std::pair<int,int> start_coordinates, int orig_value, int dest_value);

  bool InCube(const geom::Vec3& cube_center,const geom::Vec3& u0,const geom::Vec3& u1, const geom::Vec3& u2);

  std::vector<std::vector<std::vector<int> > > grid_;
  geom::Vec3 origin_;
  geom::Vec3 extent_;
  Real bin_size_;
  ost::mol::EntityView view_;
  ost::mol::EntityView solvated_view_;
  ost::mol::SurfaceHandle surf_;
  ost::mol::SurfaceVertexIDList solvated_vertices_;
  ost::mol::SurfaceHandle solvated_surf_;
  bool is_flooded_;
  int xbins_;
  int ybins_;
  int zbins_;
};


geom::Vec3 GetTiltAxis(geom::Vec3& axis);

struct EnergyF {
  EnergyF(std::vector<geom::Vec3>& p, std::vector<Real>& t_e, Real l, Real o):
          positions(p),transfer_energies(t_e), axis(geom::Vec3(0.0,0.0,1.0)),
          tilt_axis(geom::Vec3(1.0,0.0,0.0)), lambda(l), offset(o) { }

  typedef Eigen::Matrix<Real, 4, 1> XMatrixType;
  typedef Eigen::Matrix<Real, 1, 1> FMatrixType;

  Eigen::Matrix<Real,1,1> operator()(const Eigen::Matrix<Real, 4, 1>& x) const;

  geom::Vec3 axis;
  geom::Vec3 tilt_axis;
  Real lambda;
  std::vector<geom::Vec3> positions;
  std::vector<Real> transfer_energies;

  //define an offset parameter...
  //The levenberg-marquardt algorithm is designed to minimize an error function,
  //aims to converge towards zero. The offset parameter gets added at the end of
  //the energy calculation to get a positive result => forces algorithm to minimize

  Real offset;
};

struct EnergyDF {

   EnergyDF(const EnergyF& f): function(f),d_tilt(0.02),d_angle(0.02),
                         d_width(0.4),d_pos(0.4) { }

   Eigen::Matrix<Real,1,4> operator()(const Eigen::Matrix<Real, 4, 1>& x) const;

   EnergyF function;
   Real d_tilt;
   Real d_angle;
   Real d_width;
   Real d_pos;
};

ost::mol::EntityHandle FillMembraneDummies(const geom::AlignedCuboid& cuboid, const ost::mol::SurfaceHandle& surf, Real solvation_grid_bin_size, Real density);



FindMemParam MinimizeAlongZ(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies,geom::Vec3& axis, Real lambda);

geom::Vec3 RotateAroundAxis(geom::Vec3 point, geom::Vec3 axis, Real angle);

std::pair<std::pair<Real,Real>, Real> ScanAxis(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies, geom::Vec3& axis);

FindMemParam FindMembrane(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, std::vector<Real>& asa);

ost::mol::SurfaceHandle TransformSurface(ost::mol::SurfaceHandle& s, geom::Transform& t);

