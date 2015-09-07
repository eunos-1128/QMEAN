from django import forms
from django.forms import widgets
from django.core.exceptions import ValidationError
from django.conf import settings
import os
from ost.seq import CreateSequenceList, CreateSequence
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
		print "We need to check sequences and structures here"
		print structures, sequence
	

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
