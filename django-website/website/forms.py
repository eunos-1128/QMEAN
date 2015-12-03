from django import forms
from django.forms import widgets
from django.core.exceptions import ValidationError
from django.conf import settings
import os
try:
	from ost.seq import SequenceFromChain, CreateSequenceList, CreateSequence
	from ost.seq.alg import AlignToSEQRES
	from ost.io import SequenceListFromString, LoadPDB
except Exception, e:
	print e
	print 'OST not in path'
	
class UploadForm(forms.Form):
	structureUploaded = forms.CharField()
	project_name = forms.CharField(label='Project Name (Optional)', max_length=100, required=False)
	email = forms.EmailField(label='Email (Optional)',   max_length=100, required=False)
	sequence = forms.CharField(label='Sequence (SEQRES)',
							   widget=forms.Textarea(attrs={'rows':8}),required=False)
	qmeandisco = forms.BooleanField(label='QMEANDisCo',required=False)
	qmeanbrane = forms.BooleanField(label='QMEANBrane',required=False)

	def __init__(self, request, *args, **kwargs):
		self.request = request
		super(UploadForm, self).__init__(*args, **kwargs)

		for field_name, field in self.fields.items():
			if type(field.widget) in [widgets.TextInput,widgets.Select,widgets.EmailInput]:
				field.widget.attrs['class'] = 'form-control'


	def clean(self):
		cleaned_data = super(UploadForm, self).clean()
		structures = cleaned_data.get('structureUploaded')
		seq_list = cleaned_data.get('sequence')

		validated_structures = []

		#case 1: If the user did not supply seqres, then we must
		# go through every model and extract the seqres.
		if seq_list is None or len(seq_list)==0:
			for structure in structures:
				seqlist = self.extract_seqres(structure)
				for seq in seqlist:
					structure['seqres'].append({'name':seq.GetName(),'sequence':seq.string})
				validated_structures.append(structure)

			if len(validated_structures)==0:
				raise ValidationError("No valid structures were found in uploaded data!")

		#case 2: One sequence was uploaded. We must check that 
		# every chain in every model passes the AlignToSEQRES check. 
		elif len(seq_list)==1:
			for structure in structures:
				model = structure['model'].Select("peptide=true")
				for ch in model.chains:
					seqresIdentical = True
					try:
						aln = AlignToSEQRES(ch,seq_list[0].GetString(), validate=False)
						structure['seqres'].append({'name':ch.GetName(),'sequence':seq_list[0].GetString()})
					except Exception, e:
						seqresIdentical = False
				
				if seqresIdentical:
					validated_structures.append(structure)

			if len(validated_structures)==0:
				if len(structures)==1:
					if len(structures[0]['model'].chains)==1:
						raise ValidationError("Could not align the given sequence with the model!")
					else:
						raise ValidationError("Could not align the given sequence with every chain in the model!")
				else:
					raise ValidationError("Could not successfully align the given sequence with any of the models!")


		#case 3: Multiple sequences were uploaded. We have to check
		#every chain has a matching sequence based on chain/sequence name, 
		else:
			for structure in structures:
				model = structure['model'].Select("peptide=true")
				found_all_chains = True
				seqresIdentical = True
				for ch in model.chains:
					found_sequence = False
					for seq_handle in seq_list:
						if ch.GetName() == seq_handle.GetName().strip():
							found_sequence = True
							try:
								aln = AlignToSEQRES(ch,seq_handle.GetString(),validate=False)
								structure['seqres'].append({'name':ch.GetName(),'sequence':seq_handle.GetString()})
								break
							except Exception, e:
								seqresIdentical = False
					if not found_sequence:
						found_all_chains = False

				if found_all_chains and seqresIdentical:
					validated_structures.append(structure)

			if len(validated_structures)==0:
				raise ValidationError("Could not align every given sequence by chain name to the model chain name!")


		self.cleaned_data['structureUploaded'] = validated_structures



	def clean_project_name(self):
		data = self.cleaned_data['project_name']
		data = ' '.join(data.splitlines())
		return data

	def clean_structureUploaded(self):
		data = self.cleaned_data['structureUploaded']
		data = data.split(',')
		models = []
		for d in data:			
			tmppath = str(os.path.join(settings.TMP_DIR,d))
			if not os.path.exists(tmppath):
				continue
			original_name = self.request.session.get('uploaded_'+d, False)
			if not original_name:
				continue
			try:
				model = LoadPDB(tmppath)
				models.append({'tmpname':d,'tmppath':tmppath,'model':model,'original_name':original_name,'seqres':[]})
			except Exception, e:
				print e
				try:
					os.remove(os.path.join(settings.TMP_DIR,tmp_path))
				except e2:
					print e2


		if len(models)<1:
			raise ValidationError('Uploaded file no longer in session!!')
                
		return models

	def clean_sequence(self):
		#If this method is called, the user entered something in the
		#seqres box. First, can we parse this information?
		data = str(self.cleaned_data['sequence']).strip()
		seq_list = CreateSequenceList()

		if len(data)>0:
			try:
				seq_list = SequenceListFromString(data,'fasta')
			except:
				try:
					seq_list = SequenceListFromString(data,'clustal')
				except:
					try:
						seq_list.AddSequence(CreateSequence("unnamed",data))
					except Exception, e:
						pass
			for s in seq_list:
				s.SetName(s.GetName().strip())

		return seq_list

	
	def extract_seqres(self,structure):
		seq_list = CreateSequenceList()
		try:
			model = structure['model'].Select("peptide=true")
			for chain in model.GetChainList():
				seqres = ''.join([r.one_letter_code for r in chain.residues]) 
				seq_list.AddSequence(CreateSequence(chain.GetName(), seqres))
		except Exception, e:
			print 'Problem extracting seqres',e
			return None

		return seq_list
