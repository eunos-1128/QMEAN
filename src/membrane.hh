#include <ost/mol/mol.hh>

#include <ost/geom/geom.hh>
#include <vector>
#include <numeric>
#include <limits>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <qmean/module_config.hh>
#include <ost/mol/surface_handle.hh>
#include <ost/img/image_handle.hh>
#include <math.h>
#include <ost/img/image_factory.hh>
#include <ost/img/point.hh>
#include <ost/img/size.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/geom/transform.hh>

#include <ost/img/alg/levenberg_marquardt.h>


#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/SVD>
#include <Eigen/LU>


class SolvationGrid;
typedef boost::shared_ptr<SolvationGrid> SolvationGridPtr;

struct FindMemParam{
  FindMemParam() { }

  geom::Vec3 GetMembraneAxis(){
    geom::Vec3 transformed_x = transform.Apply(geom::Vec3(1.0,0.0,0.0));
    geom::Vec3 transformed_y = transform.Apply(geom::Vec3(0.0,1.0,0.0));
    return geom::Normalize(geom::Cross(transformed_x,transformed_y));
    //return transform.Apply(geom::Vec3(0.0,0.0,1.0));
  }
  Real tilt;
  Real angle;
  Real width;
  Real pos;
  geom::Transform transform;
  Real energy;
};

class DLLEXPORT_QMEAN SolvationGrid {
public:
  SolvationGrid() {};

  SolvationGrid(geom::AlignedCuboid cub, Real bin_size);

  void AddView(const ost::mol::EntityView& view) { view_ = view; };

  void AddSurface(const ost::mol::SurfaceHandle& surface);

  void Flood();

  bool IsFilled(geom::Vec3& pos);

  geom::Vec3 GetOrigin() { return origin_; }

  ost::img::ImageHandle AsIMG();

  void CalculateSolvatedSurface();

  void CalculateSolvatedView();

  ost::mol::SurfaceHandle GetSolvatedSurface();

  ost::mol::EntityView GetSolvatedView();


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
  ost::mol::SurfaceHandle solvated_surf_;
  bool is_flooded_;
  int xbins_;
  int ybins_;
  int zbins_;
};


geom::Vec3 GetTiltAxis(geom::Vec3& axis);

struct EnergyF {
  EnergyF(std::vector<geom::Vec3>& p, std::vector<Real>& t_e, Real l):
          positions(p),transfer_energies(t_e), axis(geom::Vec3(0.0,0.0,1.0)),
          tilt_axis(geom::Vec3(1.0,0.0,0.0)), lambda(l) { }

  typedef Eigen::Matrix<Real, 4, 1> XMatrixType;
  typedef Eigen::Matrix<Real, 1, 1> FMatrixType;

  Eigen::Matrix<Real,1,1> operator()(const Eigen::Matrix<Real, 4, 1>& x) const;

  geom::Vec3 axis;
  geom::Vec3 tilt_axis;
  Real lambda;
  std::vector<geom::Vec3> positions;
  std::vector<Real> transfer_energies;
};

struct EnergyDF {

   EnergyDF(const EnergyF& f): function(f),d_tilt(0.05),d_angle(0.05),
                         d_width(0.1),d_pos(0.1) { }

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

std::pair<std::pair<int,int>, Real> ScanAxis(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies, geom::Vec3& axis, Real lambda);

FindMemParam FindMembrane(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, std::vector<Real>& asa);

ost::mol::SurfaceHandle TransformSurface(ost::mol::SurfaceHandle& s, geom::Transform& t);
