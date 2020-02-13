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

