#ifndef SPHERICAL_SMOOTHER_HH
#define SPHERICAL_SMOOTHER_HH

#include <qmean/module_config.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/spatial_organizer.hh>
#include <vector>
#include <limits>
#include <math.h>
#include <ost/io/io_exception.hh>

namespace qmean{

class SphericalSmoother;
typedef boost::shared_ptr<SphericalSmoother> SphericalSmootherPtr;


struct DLLEXPORT_QMEAN SpatialOrganizerItemSmoother{

  SpatialOrganizerItemSmoother() { }
  SpatialOrganizerItemSmoother(geom::Vec3 p, int i): pos(p), index(i) { }

  geom::Vec3 pos;
  int index;

};

class DLLEXPORT_QMEAN SphericalSmoother{

public:

  SphericalSmoother() { }

  SphericalSmoother(std::vector<geom::Vec3>& positions, Real std);

  std::vector<Real> Smooth(std::vector<Real>& values);


private:

  Real Evaluate(std::vector<std::pair<Real,Real> >& pair_vec);

  std::vector<std::vector<std::pair<Real,Real> > > visitor_pattern_;
  std::vector<std::vector<int> > value_distributor_;
};

}
#endif

