from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.conf import settings
from .forms import UploadForm
import io, os, json, tempfile, tarfile, random, re, traceback, codecs
from ost.seq import CreateSequenceList, CreateSequence
from ost.io import LoadPDB, SequenceListFromString 
from run_qmean_task import RunQMEAN

def home(request):
	project_creation_error = False
	if request.method == 'POST':
		form = UploadForm(request, request.POST)
		if form.is_valid():
			project_id = create_project(request, form)
			if project_id:
				split_path = re.findall('..',project_id)
				split_path[-1] += '.qm'
				project_path = os.path.join(settings.PROJECT_DIR,os.path.sep.join(split_path))
				print "project path: ",project_path	
				qm_runner = RunQMEAN(project_path)
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

				data = {
						'meta':{},
						'options':{},
						'sequences':[],
						'models':[]
				}
				data['options']['qmeandisco'] = (True if form.cleaned_data['QMEANDisCo'] else False)
				data['options']['qmeanbrane'] = (True if form.cleaned_data['QMEANBrane'] else False)
				if form.cleaned_data['email']:
					data['meta']['email'] = form.cleaned_data['email']
				if form.cleaned_data['project_name']:
					data['meta']['project_name'] = form.cleaned_data['project_name']

				data['models'] = []
				for i,model in enumerate(form.cleaned_data['structureUploaded']):
					model_id = 'model_%03d.pdb' % (i+1)
					os.rename(model['tmppath'], os.path.join(input_dir,model_id) )
					data['models'].append({model_id:request.session['uploaded_'+model['tmpname']]})
				
				for seq in form.cleaned_data['sequence']:
					data['sequences'].append( {'name':seq.GetName(), 'sequence':seq.GetString() } )

				with io.open(os.path.join(input_dir,'project.json'), 'w', encoding='utf8') as json_file:
					data = json.dumps(data, ensure_ascii=False, encoding='utf8',indent=4, separators=(',', ': '))
					json_file.write(unicode(data))

				f = open(os.path.join(new_proj_dir,'status'),'w')
				f.write('INITIALISING')
				f.close()

				return project_id

			except Exception, e:
				print 'Failed to make new project directory %s' % project_path
				print traceback.print_exc()
				return None
			break


def sequence_upload(request):
	ret_obj = {}
	try:
		if request.FILES:
			uploaded = {"files":[] }
			seq_list = CreateSequenceList()
			for key in request.FILES:
				f = request.FILES[key]
				content = str(f.read())
				try:
                        		seq_list = SequenceListFromString(content,'fasta')
                		except:
                        		try:
                                		seq_list = SequenceListFromString(content,'clustal')
                        		except:
                                		try:
                                        		seq_list.AddSequence(CreateSequence("unnamed",''.join(content.splitlines())))
                                		except Exception, e:
                                        		print e

			if len(seq_list) > 0:
				uploaded['files'].append({"name":f.name,"content":content})
			else:
				uploaded['files'].append({"name":f.name,"error":"Could not load sequence from file"})
			ret_obj = uploaded
			
	except Exception, e:
		print e
		print traceback.print_exc()
		ret_obj['error'] = 'Sorry, can only read FASTA, Clustal or raw sequence format.'
	
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
		LoadPDB(tmp_path)
		request.session['uploaded_'+basename] = original_name
	except Exception, e:
		try:
			os.remove(os.path.join(settings.TMP_DIR,tmp_path))
		except e2:
			print e2
		err = str(e)
		tmpname=""

	return {"original_name":original_name,"file_name":basename,"error":err}


def results(request, projectid):
	status = get_project_status(projectid)
	input_data = get_project_input(projectid)
	if status in ['INITIALISING','QUEUING','RUNNING']:
		return render(request, 'results_running.html',{'projectid':projectid,
														'input_data':input_data,
														'status':status})
	return render(request,
		 			'results_completed.html' if status=='COMPLETED' else 'results_failed.html',
		 			{'projectid':projectid,
		 			'input_data':input_data,
		 			'status':status})


def get_project_input(projectid):
	try:
		f = open(os.path.join(project_path(projectid),'input','project.json'))
		data = json.load(f)
		f.close()
		return data
	except Exception, e:
		print e

def get_project_status_json(request,projectid):
	return HttpResponse(json.dumps({'status':get_project_status(projectid)}), content_type='application/json') 
	

def get_project_status(projectid):
	try:
		f = open(os.path.join(project_path(projectid),'status'))
		status = f.read().strip()
		f.close()
		return 'RUNNING'
		if status in ['INITIALISING','QUEUING','RUNNING','COMPLETED','FAILED']:
			return status
	except Exception, e:
		print e

def uploaded_structure(request, projectid, modelid):
	pdbfile=open(os.path.join(project_path(projectid),'input',modelid), 'r')
	return HttpResponse(pdbfile, content_type='text/plain;')


def project_path(projectid):
	split_path = os.path.sep.join(re.findall('..',projectid))+'.qm'
	return os.path.join(settings.PROJECT_DIR, split_path)

def help(request):
	return render(request, 'help.html')

def references(request):
	return render(request, 'references.html')

def contact(request):
	return render(request, 'contact.html')


