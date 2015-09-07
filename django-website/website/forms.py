from django import forms
from django.forms import widgets
from django.core.exceptions import ValidationError
from django.conf import settings
import os
from ost.io import SequenceFromString, LoadPDB
from ost.seq import CreateSequence
import ost

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
					model = ost.io.LoadPDB(os.path.join(str(settings.TMP_DIR),str(d)))

					seq_list = self.cleaned_data['sequence']

					print "yes, i have the seq list"
					
				except Exception, e:
					print e
					
					data.remove(d)
				

		if len(data)<1:
			raise ValidationError('Uploaded file no longer in session!!')

                
		return data

	def clean_sequence(self):
		data = str(self.cleaned_data['sequence'])
		seq = ost.seq.CreateSequenceList()
		try:
			seq = ost.io.SequenceListFromString(data,'fasta')
		except Exception, e:
			try:
				seq = ost.io.SequenceListFromString(data,'clustal')
			except Exception, e:
				try:
					ost.seq.AddSequence(ost.seq.CreateSequence("unnamed",data))
				except Exception, e:
					print e


		for s in seq:
			s.SetName(s.GetName().strip())

		if len(seq) == 0:
			raise ValidationError('Sequence input could not be read!!')
		
		return seq
