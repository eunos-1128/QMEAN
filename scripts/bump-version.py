#!/usr/bin/env python
import sys

if len(sys.argv) < 3:
  print "USAGE: python scripts/bump-version.py QMEAN_VERSION OST_VERSION"
  print "-> *_VERSION format is MAJOR.MINOR.PATCH (e.g. 1.9.1)"
  print "-> assumption is that git tags will exist for those *_VERSION"
  sys.exit(1)

# split up version number
version_string = sys.argv[1]
version = version_string.split('.')
major, minor, patch = (int(version[0]), int(version[1]), int(version[2]))
ost_version_string = sys.argv[2]

# fix CMakeLists
lines = open("CMakeLists.txt").readlines()
for i, line in enumerate(lines):
  if line.startswith("set (QMEAN_VERSION_MAJOR"):
    lines[i] = "set (QMEAN_VERSION_MAJOR %d)\n" % major
  elif line.startswith("set (QMEAN_VERSION_MINOR"):
    lines[i] = "set (QMEAN_VERSION_MINOR %d)\n" % minor
  elif line.startswith("set (QMEAN_VERSION_PATCH"):
    lines[i] = "set (QMEAN_VERSION_PATCH %d)\n" % patch
  elif line.startswith("find_package(OPENSTRUCTURE "):
    lines[i] = "find_package(OPENSTRUCTURE %s REQUIRED\n" % ost_version_string
open("CMakeLists.txt", "w").writelines(lines)

# fix CHANGELOG
lines = open("CHANGELOG.txt").readlines()
for i, line in enumerate(lines):
  if line.startswith("Changes in Release") and "X" in line.upper():
    lines[i] = "Changes in Release %s\n" % version_string
open("CHANGELOG.txt", "w").writelines(lines)

# fix doc config
lines = open("doc/source/conf.py").readlines()
for i, line in enumerate(lines):
  if line.startswith("version = "):
    lines[i] = "version = '%d.%d'\n" % (major, minor)
  elif line.startswith("release = "):
    lines[i] = "release = '%s'\n" % version_string
open("doc/source/conf.py", "w").writelines(lines)
