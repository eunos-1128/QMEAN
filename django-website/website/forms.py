from django import forms
from django.forms import widgets
from django.core.exceptions import ValidationError
from django.conf import settings
import os
try:
	from ost.seq import CreateSequenceList, CreateSequence
	from ost.seq.alg import AlignToSEQRES
	from ost.io import SequenceListFromString, LoadPDB
except:
	print 'OST not in path'
	
class UploadForm(forms.Form):
	structureUploaded = forms.CharField()
	project_name = forms.CharField(label='Project Name (Optional)', max_length=100, required=False)
	email = forms.EmailField(label='Email (Optional)',   max_length=100, required=False)
	sequence = forms.CharField(label='Sequence (SEQRES)',
							   widget=forms.Textarea)
	QMEANDisCo = forms.BooleanField(label='QMEANDisCo',required=False)
	QMEANBrane = forms.BooleanField(label='QMEANBrane',required=False)

	def __init__(self, request, *args, **kwargs):
		self.request = request
		super(UploadForm, self).__init__(*args, **kwargs)

		for field_name, field in self.fields.items():
			if type(field.widget) in [widgets.TextInput,widgets.Select,widgets.EmailInput]:
				field.widget.attrs['class'] = 'form-control'


	def clean(self):
		cleaned_data = super(UploadForm, self).clean()
		structures = cleaned_data.get('structureUploaded')
		sequence = cleaned_data.get('sequence')

		#case 1: There is only one sequence... either there are single
		#        chain models or homo-oligomers allowed
		if len(sequence) == 1:
			for s in structures:
				model = s['model'].Select("peptide=true")
				for ch in model.chains:
					try:
						aln = AlignToSEQRES(ch,sequence[0].GetString())
					except Exception, e:
						raise ValidationError("Could not align structural data to provided sequence!")

		#case 2: There is more than one sequence... For every chain we try to find a matching
		#        sequence based on chain/sequence name
		else:
			for s in structures:
				model = s['model'].Select("peptide=true")
				for ch in model.chains:
					found_sequence = False
					for seq_handle in sequence:
						#WHY IS THIS TO UPPER, LOWERCASE MAY BE IMPORTANT??
						if ch.GetName().upper() == seq_handle.GetName().upper().strip():
							found_sequence = True
							try:
								aln = AlignToSEQRES(ch,seq_handle.GetString())
								break
							except Exception, e:
								raise ValidationError("Could not align structural data to provided sequence!")
					if not found_sequence:
						raise ValidationError("Could not find an appropriate sequence for every chain!")
	


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
			if not self.request.session.get('uploaded_'+d, False):
				continue
			try:
				model = LoadPDB(tmppath)
				models.append({'tmpname':d,'tmppath':tmppath,'model':model})
			except Exception, e:
				print e

		if len(models)<1:
			raise ValidationError('Uploaded file no longer in session!!')
                
		return models

	def clean_sequence(self):
		data = str(self.cleaned_data['sequence'])
		seq_list = CreateSequenceList()
		try:
			seq_list = SequenceListFromString(data,'fasta')
		except:
			try:
				seq_list = SequenceListFromString(data,'clustal')
			except:
				try:
					seq_list.AddSequence(CreateSequence("unnamed",data))
				except:
					pass

		for s in seq_list:
			s.SetName(s.GetName().strip())

		if len(seq_list) == 0:
			raise ValidationError("Sequence input results in 0 read sequences. Please use proper FASTA format an provide at least one sequence!")

		
		return seq_list
