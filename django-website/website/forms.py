from django import forms
from django.forms import widgets
from django.core.exceptions import ValidationError
from django.conf import settings
import os
from ost.seq import CreateSequenceList, CreateSequence
from ost.seq.alg import AlignToSEQRES
from ost.io import SequenceListFromString, LoadPDB

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

		#case 1: There is no sequence... Raise an error!
		if len(sequence) == 0:
			raise ValidationError("Sequence input results in 0 read sequences. Please use proper FASTA format an provide at least one sequence!")

		#case 2: There is only one sequence... either there are single
		#        chain models or homo-oligomers allowed
		elif len(sequence) == 1:
			for s in structures:
				model = LoadPDB(os.path.join(settings.TMP_DIR,str(s))).Select("peptide=true")
				for ch in model.chains:
					try:
						aln = AlignToSEQRES(ch,sequence[0].GetString())
					except Exception, e:
						raise ValidationError("Could not align structural data to provided sequence!")

		#case 3: There is more than one sequence... For every chain we try to find a matching
		#        sequence based on chain/sequence name
		else:
			for s in structures:
				model =  LoadPDB(os.path.join(settings.TMP_DIR,str(s))).Select("peptide=true")
				for ch in model.chains:
					found_sequence = False
					for seq_handle in sequence:
						print ch.GetName(),seq_handle.GetName()
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
		for d in data[:]:
			removed = False			
			if not os.path.exists(os.path.join(str(settings.TMP_DIR),str(d))):
				data.remove(d)
				removed = True
			if not self.request.session.get('uploaded_'+d, False):
				data.remove(d)
				removed = True
			if not removed:
				try:
					model = LoadPDB(os.path.join(str(settings.TMP_DIR),str(d)))
				except Exception, e:
					data.remove(d)

		if len(data)<1:
			raise ValidationError('Uploaded file no longer in session!!')
                
		return data

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
			raise ValidationError('Sequence input could not be read!!')
		
		return seq_list
