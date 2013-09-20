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
