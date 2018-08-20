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

#include <boost/python.hpp>
#include <ost/mol/entity_handle.hh>

void export_Torsion();
void export_Interaction();
void export_Packing();
void export_Reduced();
void export_CBeta();
void export_SSAgreement();
void export_Base();
void export_SphericalSmoother();
void export_Membrane();
void export_CBPacking();
void export_HBond();
void export_disco();
void export_clash();

using namespace boost::python;

BOOST_PYTHON_MODULE(_qmean)
{
  export_Base();	
  export_Torsion();
  export_Interaction();
  export_Packing();
  export_Reduced();
  export_CBeta();
  export_SSAgreement();
  export_SphericalSmoother();
  export_Membrane();
  export_CBPacking();
  export_HBond();
  export_disco();
  export_clash();
}

