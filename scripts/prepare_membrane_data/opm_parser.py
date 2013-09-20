"""
takes opm data as input. The datafile can be 
downloaded at http://opm.phar.umich.edu/subunits.php
the output is a dictionary with pdb_ids as keys and
ost queries as values

parameters:
  opm_data:     read in above mentioned file with open(file,'r').readlines()
                and pass it as input
  selected_ids: list of strings. only structures with pdb-id present in 
                this list will be added to output.
  only_chains:  if set to True, queries only consiting of chains will be generated,
                such as: 'cname=A or cname=B'
  save:         possibility to pass a filename. If filename is set, the output will be
                written in this file with format: pdbid;query
                =>one line per query
output:
  dictionary with pdb_ids as keys and queries as values


"""

def ParseOPM(opm_data, selected_ids=None, only_chains=False, save=False):

  import re

  queries=dict()
  old_id=''
  query_contributions=list()
  rnum_regex=re.compile(r'\d+\s*\(\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)\s*\)')
  #add additional line, otherwise the last query wouldn't be added
  opm_data.append('end of file')

  for line in opm_data:
    current_id=line.split()[0]
    if current_id!=old_id:
      if len(query_contributions)>0:
        if selected_ids:
          if old_id in selected_ids:
            queries[old_id]=' or '.join(query_contributions)
        else:
          queries[old_id]=' or '.join(query_contributions)
      old_id=current_id
      query_contributions=list()

    cname=line.split()[1]

    if only_chains:
      query_contributions.append('cname='+cname)
      continue

    rnum_groups=re.findall(rnum_regex, line)
    temp=list()

    for item in rnum_groups:
      temp.append(item[0]+':'+item[1])

    query_contributions.append('(cname='+cname+' and rnum='+','.join(temp)+')')

  if save:
    outfile=open(save,'w')
    for k,v in queries.iteritems():
      outfile.write(k+';'+v)
      outfile.write('\n')
    outfile.close()
  return queries 
   