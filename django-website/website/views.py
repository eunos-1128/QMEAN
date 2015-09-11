from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.conf import settings
from django.core.servers.basehttp import FileWrapper
from django.contrib.staticfiles.finders import find
from .forms import UploadForm
import io, os, json, tempfile, tarfile, random, re, traceback, codecs, zipfile, glob
from datetime import datetime
try:
	from ost.seq import CreateSequenceList, CreateSequence
	from ost.io import LoadPDB, SequenceListFromString 
	from run_qmean_task import RunQMEAN
except:
	print 'OST not in path'

def home(request):
	project_creation_error = False
	if request.method == 'POST':
		form = UploadForm(request, request.POST)
		if form.is_valid():
			project_id = create_project(request, form)
			if project_id:
				qm_runner = RunQMEAN(project_path(project_id))
				return HttpResponseRedirect(reverse('results',args=[project_id]))
			else:
				project_creation_error = 'Failed to create a new project'
	else:
		form = UploadForm(request)

	return render(request, 'home.html', {'form':form, 'project_creation_error':project_creation_error})

def create_project(request, form):

	for attempt in range(10):
		project_id = ''.join([random.choice('abcdefghjkmnpqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ23456789') for x in range(6)]) 
		new_project_path = project_path(project_id)

		if not os.path.exists(new_project_path):
			try:
				os.makedirs(new_project_path)
				input_dir = os.path.join(new_project_path,'input')
				os.makedirs(input_dir)
				os.makedirs(os.path.join(new_project_path,'output'))

				data = {
						'meta':{'created':str(datetime.now())},
						'options':{},
						'sequences':[],
						'models':[]
				}
				data['options']['qmeandisco'] = (True if form.cleaned_data['QMEANDisCo'] else False)
				if form.cleaned_data['email']:
					data['meta']['email'] = form.cleaned_data['email']
				if form.cleaned_data['project_name']:
					data['meta']['project_name'] = form.cleaned_data['project_name']

				data['models'] = []
				for i,model in enumerate(form.cleaned_data['structureUploaded']):
					modelid = 'model_%03d' % (i+1)
					os.rename(model['tmppath'], os.path.join(input_dir,modelid+'.pdb') )
					data['models'].append({'modelid':modelid,'name':request.session['uploaded_'+model['tmpname']]})
				
				for seq in form.cleaned_data['sequence']:
					data['sequences'].append( {'name':seq.GetName(), 'sequence':seq.GetString() } )

				with io.open(os.path.join(input_dir,'project.json'), 'w', encoding='utf8') as json_file:
					data = json.dumps(data, ensure_ascii=False, encoding='utf8',indent=4, separators=(',', ': '))
					json_file.write(unicode(data))

				f = open(os.path.join(new_project_path,'status'),'w')
				f.write('INITIALISING')
				f.close()

				return project_id

			except Exception, e:
				print 'Failed to make new project directory %s' % new_project_path
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

	return {"original_name":original_name,"file_name":basename,"error":err}


def results(request, projectid):
	status = get_project_status(projectid)
	input_data = get_project_input(projectid)
	if 'created' in input_data['meta'].keys():
		input_data['meta']['created'] = datetime.strptime(input_data['meta']['created'],"%Y-%m-%d %H:%M:%S.%f")
	if status in ['INITIALISING','QUEUEING','RUNNING']:
		return render(request, 'results_running.html',{'projectid':projectid,
														'input_data':input_data,
														'status':status})

	for model in input_data['models']:
		try:
			global_scores = {}
			mdl_dir = os.path.join(project_path(projectid),'output',model['modelid'])
			f = open(os.path.join(mdl_dir,'global_scores.txt'))
			lines = f.readlines()
			f.close()
			for line in lines:
				if line.startswith('#') or line.startswith('name'):
					continue
				scores = line.split()
				global_scores[scores[0]] = {'norm':scores[1],'zscore':scores[2]}
			model['global_scores'] = global_scores

			model['local_quality_plots'] = []
			for f in glob.iglob(os.path.join(mdl_dir,'plots','local_quality_estimate_*')):
				if os.path.isfile(f):
					f = f.split(os.path.sep)[-1]
					m = re.search('_([\w])\.png',f)
					if m:
						model['local_quality_plots'].append(f)



		except Exception, e:
			print e

	input_data['models'].sort( key=lambda x: x['global_scores']['qmean4'] )

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
		status = f.read().strip().split()[-1]
		f.close()
		if status in ['INITIALISING','QUEUEING','RUNNING','COMPLETED','FAILED']:
			return status
	except Exception, e:
		print e

def model_structure(request, projectid, modelid):
	pdbfile=open(os.path.join(project_path(projectid),'output',modelid,'model.pdb'), 'r')
	return HttpResponse(pdbfile, content_type='text/plain;')

def uploaded_structure(request, projectid, modelid):
	pdbfile=open(os.path.join(project_path(projectid),'input',modelid+'.pdb'), 'r')
	return HttpResponse(pdbfile, content_type='text/plain;')

def archive(request, projectid, modelid=None):

	if get_project_status(projectid)!='COMPLETED' and get_project_status(projectid)!='FAILED':
		return HttpResponse('Please wait until template search has completed!<br><br><a href="#" onclick="window.close()">Close this window</a>')

	ppath = project_path(projectid)

	input_data = get_project_input(projectid)	
	if 'project_name' in input_data['meta'].keys():
		archivename = input_data['meta']['project_name'] 
	else:
		archivename = 'QMEAN_Project'

	archivename += '_'+ input_data['meta']['created'].split()[0] 
	safearchivename = []
	for c in archivename:
		if re.match('[\da-zA-Z _|]',c):
			safearchivename.append(c)
	archivename = ''.join(safearchivename)
	archivename = archivename.replace(' ','_')

	if modelid is None:
		archive_path = os.path.join(ppath,'archive.zip')
	else:
		archive_path = os.path.join(ppath,'output',modelid,'archive.zip')

	# if archive file has been created, just return here.    
	if os.path.isfile(archive_path):
		archive_file = open(archive_path,'r')
		wrapper = FileWrapper(archive_file)
		response = HttpResponse(wrapper, content_type='application/zip')
		response['Content-Disposition'] = 'attachment; filename=%s.zip' %(archivename)
		response['Content-Length'] = archive_file.tell()
		archive_file.seek(0)
		return response


	archive_file = open(archive_path,'w')
	archive = zipfile.ZipFile(archive_file, 'w', zipfile.ZIP_DEFLATED) 

	archive.write( os.path.join(project_path(projectid), 'input', 'project.json'), 'project.json')   
	    
	for m in input_data['models']:
		if modelid is not None and modelid!=m['modelid']:
			continue
		mdl_dir = os.path.join(project_path(projectid), 'output', m['modelid']) 

		safemodelname = []
		for c in str(m['name']).split('.pdb')[0]:
			if re.match('[\da-zA-Z _|]',c):
				safemodelname.append(c)
		safemodelname = ''.join(safemodelname)

		archive.write( os.path.join(project_path(projectid), 'input', m['modelid']+'.pdb'),
						 safemodelname+os.path.sep+m['name'])


		for f in ['model.pdb','global_scores.txt','local_scores.txt']:
			if os.path.exists(os.path.join(mdl_dir,f)):
				archive.write( os.path.join(mdl_dir,f), safemodelname+os.path.sep+f)
  
		for f in glob.iglob(os.path.join(mdl_dir,'plots','*')):
			if os.path.isfile(f):
				archive.write( f, safemodelname+"/plots/%s" % os.path.basename(f))


	archive.close()
	archive_file.close()

	archive_file = open(archive_path,'r')
	wrapper = FileWrapper(archive_file)
	response = HttpResponse(wrapper, content_type='application/zip')
	response['Content-Disposition'] = 'attachment; filename=%s.zip' %(archivename)
	response['Content-Length'] = archive_file.tell()
	archive_file.seek(0)
	return response


def qmean_pix(request, projectid, modelid, name):

	img_filename = os.path.join(project_path(projectid),
								'output',
								modelid,
								'plots',
								name)
	if not os.path.exists(img_filename):
		img_filename = find('img/failed.jpg')

	image_data = open(img_filename, "rb").read()  
	return HttpResponse(image_data, content_type="image/png")



def project_path(projectid):
	split_path = os.path.sep.join(re.findall('..',projectid))+'.qm'
	return os.path.join(settings.PROJECT_DIR, split_path)

def help(request):
	return render(request, 'help.html')

def references(request):
	return render(request, 'references.html')

def contact(request):
	return render(request, 'contact.html')


