#include <ost/mol/mol.hh>

#include <ost/geom/geom.hh>
#include <vector>
#include <numeric>
#include <limits>
#include <qmean/module_config.hh>
#include <ost/mol/surface_handle.hh>
#include <ost/img/image_handle.hh>
#include <math.h>
#include <ost/img/image_factory.hh>
#include <ost/img/point.hh>
#include <ost/img/size.hh>
#include <ost/mol/atom_handle.hh>



class SolvationGrid;
typedef boost::shared_ptr<SolvationGrid> SolvationGridPtr;

struct FindMemParam{
  FindMemParam() { }
  FindMemParam(Real t, Real a, Real w, Real p, geom::Vec3 ax):
              tilt(t), angle(a), width(w), pos(p),axis(ax) { }
  Real tilt;
  Real angle;
  Real width;
  Real pos;
  geom::Vec3 axis;
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

//struct EnergyEvaluator {
//  float operator()(Eigen::Mat3& x) {
//
//  }
//};

ost::mol::EntityHandle FillMembraneDummies(const geom::AlignedCuboid& cuboid, const ost::mol::SurfaceHandle& surf, Real solvation_grid_bin_size, Real density);

void MinimizeAlongAxis(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies,geom::Vec3& axis);

geom::Vec3 RotateAroundAxis(geom::Vec3 point, geom::Vec3 axis, Real angle);

std::pair<std::pair<int,int>, Real> ScanAxis(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies, geom::Vec3& axis);

void FindMembrane(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, std::vector<Real>& asa);


