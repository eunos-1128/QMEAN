#include <ost/mol/mol.hh>

#include <ost/geom/geom.hh>
#include <vector>
#include <qmean/module_config.hh>
#include <ost/mol/surface_handle.hh>
#include <ost/img/image_handle.hh>
#include <math.h>
#include <ost/img/image_factory.hh>
#include <ost/img/point.hh>
#include <ost/img/size.hh>



class SolvationGrid;
typedef boost::shared_ptr<SolvationGrid> SolvationGridPtr;

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


ost::mol::EntityHandle FillMembraneDummies(const geom::AlignedCuboid& cuboid, const ost::mol::SurfaceHandle& surf, Real solvation_grid_bin_size, Real density);


