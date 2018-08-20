# QMEAN:
# ----------
# 
# Copyright (C) 2018 University of Basel and the Authors.
# Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""
Author: Marco Biasini
"""

class ModelSet:
  """
  Generator for PISCES model set. Iterates over all models listed in a PISCES
  file.
  
  Usage:
  for pdb_id, chain in ModelSet('path_to_pisces_file'):
    print pdb_id, chain
  """
  def __init__(self, pisces_list, skip_header=True, emit_chain=True):
    self.pisces_list_=pisces_list
    self.skip_header_=skip_header
    self.emit_chain_=emit_chain
  
  def __iter__(self):
    pisces=open(self.pisces_list_)
    if self.skip_header_:
      pisces.readline()
    for line in pisces:
      if line=='':
        continue
      pdb_id_and_chain=(line.split(' ')[0]).strip()
      if self.emit_chain_:
        if len(pdb_id_and_chain) == 5:
          yield pdb_id_and_chain[:4], pdb_id_and_chain[4:]
        else:
          yield pdb_id_and_chain[:4], ''
      else:
        yield pdb_id_and_chain[:4]
    pisces.close()
