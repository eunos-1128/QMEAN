from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.conf import settings
from .forms import UploadForm
import os, json, tempfile, tarfile, random, re, traceback, codecs

def home(request):
	project_creation_error = False
	if request.method == 'POST':
		form = UploadForm(request, request.POST)
		if form.is_valid():
			project_id = create_project(request, form)
			if project_id:
				return HttpResponseRedirect(reverse('results',args=[project_id]))
			else:
				project_creation_error = 'Failed to create a new project'
	else:
		form = UploadForm(request)

	return render(request, 'home.html', {'form':form, 'project_creation_error':project_creation_error})

def create_project(request, form):

	for attempt in range(10):
		project_id = ''.join([random.choice('abcdefghjkmnpqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ23456789') for x in range(6)]) 
		split_path = re.findall('..',project_id)
		split_path[-1] += '.qm'
		project_path = os.path.join(settings.PROJECT_DIR,os.path.sep.join(split_path))

		if not os.path.exists(project_path):
			try:
				new_proj_dir = settings.PROJECT_DIR
				for p in split_path:
					new_proj_dir = os.path.join(new_proj_dir,p)
					os.makedirs(new_proj_dir)
				input_dir = os.path.join(new_proj_dir,'input')
				os.makedirs(input_dir)
				os.makedirs(os.path.join(new_proj_dir,'output'))

				f = codecs.open(os.path.join(input_dir,'meta.txt'),'w','utf-8')
				f.write('QMEANDISCO\t%s\n' % ('True' if form.cleaned_data['QMEANDisCo'] else 'False'))
				f.write('QMEANBRANE\t%s\n' % ('True' if form.cleaned_data['QMEANBrane'] else 'False'))
				if form.cleaned_data['email']:
					f.write('EMAIL\t%s\n' % form.cleaned_data['email'])
				if form.cleaned_data['project_name']:
					f.write('PROJECT_NAME\t%s\n' % form.cleaned_data['project_name'])

				for i,tmpname in enumerate(form.cleaned_data['structureUploaded']):
					model_id = 'model_%03d.pdb' % (i+1)
					os.rename(os.path.join(settings.TMP_DIR,tmpname), os.path.join(input_dir,model_id) )
					f.write('%s\t%s\n' %(model_id,request.session['uploaded_'+tmpname]))
				f.close()
				if form.cleaned_data['sequence']:
					f = open(os.path.join(input_dir,'target.fasta'),'w')
					f.write(form.cleaned_data['sequence'])
					f.close()

				f = open(os.path.join(input_dir,'status'),'w')
				f.write('INITIALISING')
				f.close()

				return project_id

			except Exception, e:
				print 'Failed to make new project directory %s' % project_path
				print traceback.print_exc()
				return None
			break


def sequence_upload(request):
	if request.FILES:
		uploaded = {"files":[] }
		for key in request.FILES:
			f = request.FILES[key]
			uploaded['files'].append({"name":f.name,
									  "size":f.size,
									  "content":f.read()})
		ret_obj = uploaded
	return HttpResponse(json.dumps(ret_obj), content_type='application/json') 


def structure_upload(request):
	if request.FILES:
		uploaded = {"files":[] }
		for key in request.FILES:
			f = request.FILES[key]
			if f.size > 10485760:
				uploaded['files'].append({"name":f.name,
									  "size":f.size,
									  "error":'File too large, must be less than 10MB'})
			else:
				try:
					tmpfile = tempfile.NamedTemporaryFile(dir=settings.TMP_DIR,delete=False)
					tmp_name = tmpfile.name
					with open(tmp_name, 'wb+') as destination:
						for chunk in f.chunks():
							destination.write(chunk)
					tmpfile.close()

					if not tarfile.is_tarfile(tmp_name):
						valid_struc = is_valid_structure_file(request, f.name, tmp_name)
						uploaded['files'].append({"name":valid_struc['original_name'],
							"file_id":valid_struc['file_name'],
							"error":valid_struc['error']})
					else:
						tar = tarfile.open(tmp_name)
						for member in tar.getmembers():
							f=tar.extractfile(member)
							tmpfile = tempfile.NamedTemporaryFile(dir=settings.TMP_DIR,delete=False)
							tmp_name = tmpfile.name
							with open(tmp_name, 'wb+') as destination:
								destination.write(f.read())
							tmpfile.close()
							valid_struc = is_valid_structure_file(request, f.name, tmp_name)
							uploaded['files'].append({"name":valid_struc['original_name'],
								"file_id":valid_struc['file_name'],
								"error":valid_struc['error']})
						
						
				except Exception, e:
					print e
					uploaded['files'].append({"name":f.name,"error":e})


		ret_obj = uploaded

	return HttpResponse(json.dumps(ret_obj,indent=4, separators=(',', ': ')), content_type='application/json') 

def is_valid_structure_file(request,original_name,tmp_path):
	err = ""
	basename = os.path.basename(tmp_path)
	try:
		#This is where we will determine if the uploaded file
		#is a valid structure. If it is, we add the tmp file path
		#to the session, which ties in the users name for the file
		f = open(tmp_path,'r')
		f.close()
		request.session['uploaded_'+basename] = original_name
	except Exception, e:
		try:
			os.remove(os.path.join(settings.TMP_DIR,tmp_path))
		except e2:
			print e2
		err = str(e)
		tmpname=""

	return {"original_name":original_name,"file_name":basename,"error":err}


def results(request, jobid):
	return render(request, 'results.html')

def help(request):
	return render(request, 'help.html')

def references(request):
	return render(request, 'references.html')

def contact(request):
	return render(request, 'contact.html')

