# Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
# Biozentrum - University of Basel
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from ost import seq
from ost import mol
from qmean import *
import traceback

def AlignChainToSEQRES(chain, seqres):
  """
  aligns the chain to the seqres sequence, first by using AlignToSEQRES,
  and if that fails (usually when connectivity problems are present) 
  by aligning the sequence with Smith-Waterman using a very strict 
  substitution matrix.
  """
  try:
    return seq.alg.AlignToSEQRES(chain.Select(''), seqres, 
                                 try_resnum_first=True, validate=False)
  except Exception as e:
    print(e)

  chain_seq = ''.join([r.one_letter_code for r in chain.residues])
  aln = seq.alg.SemiGlobalAlign(seq.CreateSequence('seqres',seqres),
                                seq.CreateSequence('atom_seq',chain_seq), 
                                seq.alg.BLOSUM100)[0]

  #the seqres must not contain gaps...
  if str(aln.sequences[0]).find('-') != -1:
    raise ValueError("Could not align atomseq to seqres provided in predicted sequence feature input!")
  
  new_aln = seq.CreateAlignment(seq.CreateSequence('seqres',seqres),
                                seq.CreateSequence('atom_seq','-'*len(seqres)))

  start = aln.sequences[0].offset
  end = start + len(str(aln.sequences[0]))

  new_aln[start:end].Replace(aln[:])
  new_aln.SetSequenceOffset(1,aln.sequences[1].offset)

  return new_aln
  

class PSIPREDHandler:
  def __init__(self,psipred_data):
    try:
      self.seq = psipred_data["seq"]
      self.ss = psipred_data["ss"]
      self.conf = psipred_data["conf"]
    except:
      try:
        if not hasattr(psipred_data, 'read'):
          stream = open(psipred_data, 'r')
        else:
          stream = psipred_data
        self.seq, self.ss, self.conf = self.ParsePSIPRED(stream)
      except:
        raise ValueError("could not load psipred data!")

    if len(self.seq)!=len(self.ss):
      raise ValueError("Length of provided psipred sequence is not consistent with predicted secondary structure!")
    if len(self.seq)!=len(self.conf):
      raise ValueError("Length of provided psipred sequence is not consistent with psipred confidence!")


  def ParsePSIPRED(self,stream):
    lines=stream.readlines()
    if "# PSIPRED HFORMAT" in lines[0]:
      conf=""
      ss=""
      sequence=""
      for line in lines:
        splitted=line.strip().split() 
        if len(splitted)==2: 
          if splitted[0] == "Conf:":
            conf+=splitted[1].strip()
          if splitted[0] == "Pred:":
            ss+=splitted[1].strip()
          if splitted[0] == "AA:":
            sequence+=splitted[1].strip()
    else:
      raise RuntimeError("invalid PsiPred file!")
    return (sequence, ss, conf)

  def CutPsiPredStuff(self, chain):
    al = AlignChainToSEQRES(chain, self.seq)
    cut_ss=''
    cut_conf=''
    cut_sequence=''

    data_index=0

    for col in al:
      if col[1]!='-':
        if col[1]==col[0]:
          cut_ss+=self.ss[data_index]
          cut_conf+=str(self.conf[data_index])
          cut_sequence+=col[0]
          data_index+=1
        else:
          if col[0]!='-':
            raise ValueError("chain sequence ist not consistent with psipred sequence! ("+col[0]+"vs."+col[1]+")")
          else:
            cut_ss+='X'
            cut_conf+='X'
            cut_sequence+=col[0]
      else:
        if col[0]!='-':
          data_index+=1



    if len(chain.residues)!=len(cut_ss):
      raise RuntimeError("extracted psipred stuff has not the same length as the residues in the structure")
    if len(chain.residues)!=len(cut_conf):
      raise RuntimeError("extracted psipred stuff has not the same length as the residues in the structure")

    return (cut_ss, cut_conf, cut_sequence)

  def GetPSIPREDSS(self,chain=None):
    if chain == None:
      return self.ss
    else:
      ss, conf, seq = self.CutPsiPredStuff(chain)
      return ss
  
  def GetPSIPREDConf(self,chain=None):
    if chain == None:
      return self.conf
    else:
      ss, conf, seq = self.CutPsiPredStuff(chain)
      return conf

  @staticmethod
  def GetSSAgreement(observed_ss, psipred_ss, psipred_conf):

    if len(observed_ss)!=len(psipred_ss):
      raise ValueError("observed ss and psipred ss must have the same length")
    if len(observed_ss)!=len(psipred_conf):
      raise ValueError("observed ss and psipred conf must have the same length")

    return_list=list()
    SSCalculator=SSAgreement()

    for obs_ss, p_ss, p_c in zip(observed_ss, psipred_ss, psipred_conf):

      if p_ss == 'X' or obs_ss == 'X' or p_c=='X':
        return_list.append(float('NaN'))
        continue
      return_list.append(SSCalculator.LogOddScore(str(obs_ss),str(p_ss),int(p_c)))

    return return_list

  def GetSSAgreementFromChain(self,chain):

    psipred_ss, psipred_conf, psipred_seq = self.CutPsiPredStuff(chain)

    observed_ss=''
    for r in chain.residues:
      observed_ss+=str(r.GetSecStructure())

    return self.GetSSAgreement(observed_ss, psipred_ss, psipred_conf)


class ACCPROHandler:

  def __init__(self,accpro_data):
    try:
      self.seq=accpro_data['seq']
      self.acc=accpro_data['acc']
      try:
        self.ss=accpro_data['ss']
      except:
        self.ss=len(self.seq)*['X']
    except:
      try:
        if not hasattr(accpro_data, 'read'):
          stream = open(accpro_data, 'r')
        else:
          stream = accpro_data
        self.seq, self.ss, self.acc = self.ParseACCPRO(stream)
      except:
        print(traceback.print_exc())
        raise ValueError("could not load ACCPRO data!")

    if len(self.seq)!=len(self.ss):
      raise ValueError("length of accpro sequence is not consistent with accpro ss_prediction!")
    if len(self.seq)!=len(self.acc):
      raise ValueError("length of accpro sequence is not consistent with accpro accessiblity prediction!")

  def ParseACCPRO(self,stream):
    lines=stream.readlines()
    l=lines[1:]
    l=[i.strip() for i in l]
    #format: [seq, ss, acc]
    return l

  def CutACCPROStuff(self,chain):
    al = AlignChainToSEQRES(chain, self.seq)

    cut_ss=''
    cut_burial=''
    cut_sequence=''

    data_index=0

    for col in al:
      if col[1]!='-':
        if col[1]==col[0]:
          cut_ss+=self.ss[data_index]
          cut_burial+=self.acc[data_index]
          cut_sequence+=self.seq[data_index]
          data_index+=1
        else:
          if col[0]!='-':
            raise ValueError("chain sequence is not consistent with accpro sequence!")
          else:
            cut_ss+='X'
            cut_burial+='X'
            cut_sequence+=col[0]
      else:
        if col[0]!='-':
          data_index+=1


    return (cut_ss, cut_burial, cut_sequence)


  def GetACCPROAccessibility(self,chain=None):
    if chain == None:
      return self.acc
    else:
      accpro_ss, accpro_acc, accpro_seq = self.CutACCPROStuff(chain)
      return accpro_acc


  @staticmethod
  def GetACCAgreement(observed_acc, accpro_acc):

    return_list=list()

    if len(observed_acc)!=len(accpro_acc):
      raise ValueError("observed accessibility and accpro accessibility must have the same length!")

    for a, b in zip(observed_acc, accpro_acc):
      if a not in ['e','b'] or b not in ['e','b']:
        return_list.append(float('NaN'))
        continue
      if(a==b):
        return_list.append(1.0)
      else:
        return_list.append(0.0)

    return return_list


  def GetACCAgreementFromChain(self, chain):

    accpro_ss, accpro_acc, accpro_seq = self.CutACCPROStuff(chain)

    observed_acc = ['X'] * len(chain.residues)

    for idx, r in enumerate(chain.residues):
      try:
        if r.GetFloatProp("asaRel") < 0.25:
          observed_acc[idx] = 'b'
        else:
          observed_acc[idx] = 'e'
      except:
        # value remains 'X' => gives NaN in the end
        pass

    observed_acc = ''.join(observed_acc)

    return self.GetACCAgreement(observed_acc, accpro_acc)
