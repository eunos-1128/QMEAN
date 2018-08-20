// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

