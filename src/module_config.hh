#ifndef QMEAN_MODULE_CONFIG_HH
#define QMEAN_MODULE_CONFIG_HH

#include <ost/base.hh>

#if defined(OST_MODULE_QMEAN)
#  define DLLEXPORT_QMEAN DLLEXPORT
#else
#  define DLLEXPORT_QMEAN DLLIMPORT
#endif

#endif

