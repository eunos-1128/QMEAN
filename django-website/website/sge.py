import subprocess, shlex, threading, Queue
import qmean_server_config as config
import time, re, os, fcntl
import socket

def _IsStringLike(string):
  try:
    string+''
    return True
  except:
    return False

def _QsubTimeout(p, cmd, pipe):
  if p.poll() is None:
    try:
      p.kill()
    except:
      pass
    else:
      pipe.put('qsub did not finish in time for "%s".' % cmd)

def IsAlive(job_id):
  '''
  Check whether the job is running, queueing.

  :param job_id: Gridengine identifier.
  :type job_id: :class:`int`
  '''
  pattern='%s ' % str(job_id)
  ps = subprocess.Popen([config.QSTAT_BIN, '-u', os.environ.get("USER")],
                        shell=True,
                        stdout=subprocess.PIPE)
  for line in ps.stdout.readlines():
    if line.strip().startswith(pattern):
      return True
  return False

def Qacct(job_id):
    '''
    Get job exit status by the grid engine (not status of the job script!).
    '''
    tries = 15
    result = None
    while tries:
      ps = subprocess.Popen([config.QACCT_BIN, '-j', job_id],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      has_error = False
      for line in ps.stderr.readlines():
        if re.match(r'error: job id \d+ not found', line):
          has_error = True
          time.sleep(1)
          tries -= 1
          continue
      if not has_error:
        result = dict()
        break
    for line in ps.stdout.readlines():
      line = line.strip()
      mc = re.match(r'failed\s+(\d+)', line)
      if mc:
        result['failed'] =  int(mc.group(1))
        continue
      mc = re.match(r'exit_status\s+(\d+)', line)
      if mc:
        result['exit_status'] =  int(mc.group(1))
        continue
      mc = re.match(r'maxvmem\s+(.*)', line)
      if mc:
        result['maxvmem'] = mc.group(1)
        continue
      mc = re.match(r'maxrss\s+(.*)', line)
      if mc:
        result['maxrss'] = mc.group(1)
        continue
      mc = re.match(r'ru_wallclock\s+(\d+)', line)
      if mc:
        result['ru_wallclock'] = int(mc.group(1))
        continue
      mc = re.match(r'qname\s+(.*)', line)
      if mc:
        result['qname'] = mc.group(1)
        continue
      mc = re.match(r'hostname\s+(.*)', line)
      if mc:
        result['hostname'] = mc.group(1)
        continue
    return result

class QacctFailedError(RuntimeError):
    '''
    An exception raised when running qacct failed. That is not the job itself
    got killed by the grid engine but executing the qacct binary did not work
    out.

    :param job_id: System id of a job ran via the grid engine.
    :type job_id: :class:`int`
    '''
    def __init__(self, job_id):
        self.job_id=job_id
        RuntimeError.__init__(self, "Running qacct for job id '%s' failed." % \
                              self.job_id)

class QacctInfo(object):
    '''
    Get information via :command:`qacct` for a job. Creating an instance of
    this class will cache the results. Since instanciating only makes sense
    after a job has finished, those should not change anymore. If ran before the
    job is finished, an exception will be raised.

    Attributes of this class closely follow the names in the qacct output.

    :param job_id: System id of a job ran via the grid engine.
    :type job_id: :class:`int`
    '''
    def __init__(self, job_id):
        self.job_id = job_id
        # run qacct, distribute info
        qacct_info = Qacct(self.job_id)
        if not qacct_info:
            raise QacctFailedError(self.job_id)
        self.failed = qacct_info['failed']
        self.exit_status = qacct_info['exit_status']
        self.maxvmem = qacct_info['maxvmem']
        self.maxrss = qacct_info['maxrss']
        self.ru_wallclock = qacct_info['ru_wallclock']
        self.qname = qacct_info['qname']
        self.hostname = qacct_info['hostname']

    def GotKilled(self):
        '''
        Check wether the grid engine killed our job.
        '''
        if self.exit_status == 137: # mem
            return True
        elif self.exit_status == 138: # time
            return True
        return False

    def GotKilledOfMem(self):
        '''
        Check if a job was killed because he used more memory than allocated.
        '''
        if self.exit_status == 137:
            return True
        return False

    def GotKilledOfTime(self):
        '''
        Check if a job was killed because he used more time than allocated.
        '''
        if self.exit_status == 138:
            return True
        return False
#- job.wait() gets exceptions!

class SGEJob(object):
  '''
  Create a job for the Sun Grid Engine (SGE)
  
  :param command: The command to be executed. Either a string containing the 
      command string or a list containing the arguments.
  :param inherit_env: If set to true, the current environment is inherited for 
      the job on the cluster node
  :param queue: The queue to use for the job
  :param stdout: If set, the filename to be used for stdout of the job
  :param stderr: If set, the filename to be used for stderr of the job
  :param submit: If set to true, the job will immediately be submitted.
  :param memory: Request memory, given as tuple like (4.0, 'G').
  '''
  # info on the queues:
  # infinite.q - no runtime limitation. Access to 20% of the cpus
  # very_long.q - 1 week max runtime. Access to 60% of cpus
  # long.q - 24h max runtime. Access to 80% of the cpus
  # short.q - 6h max runtime. access to 100% of the cpus
  QUEUES = ('short.q', 'long.q', 'infinite.q', 'very_long.q',
            # The 'sync' queue is special to the productive server, its used
            # for syncing to the failover site at low priority.
            'sync')
  MT_ENVS = {'short.q': 'smp', 'long.q': 'smp', 'infinite.q': 'smp',
             'very_long.q': 'smp'}
  MT_SLOTS = {'smp': 64}
  DEFAULT_MEM = (1.0, 'G')
  def __init__(self, command, stdout=None, stderr=None, queue='long.q',
               cwd=True, name='job.sh', submit=True, inherit_env=False,
               cpu=None, hold_jid=None, memory=None):
    queue_name=queue.split('@@')[0]
    if not queue_name in SGEJob.QUEUES:
      raise ValueError('queue must be one of %s' % ', '.join(SGEJob.QUEUES))
    if cpu:
      if not SGEJob.MT_ENVS[queue_name]:
        raise ValueError('queue "%s" has no parallel environment defined' % \
                           queue_name)
      if cpu > SGEJob.MT_SLOTS[SGEJob.MT_ENVS[queue_name]]:
        raise ValueError('To many slots requested for '+
                         '"%s %s", max. %d slots, requested %d.' %
                           (queue_name, SGEJob.MT_ENVS[queue_name],
                            SGEJob.MT_SLOTS[SGEJob.MT_ENVS[queue_name]], cpu))

    self.cwd=cwd
    self.queue=queue
    self.stdout=stdout
    self.name=name
    self.stderr=stderr
    self.inherit_env=inherit_env
    self.command=command
    self.submitted=False
    self.cpu=cpu
    self.hold_jid=hold_jid
    self.memory=memory
    self._sge_info = None
    if submit:
      print "do the submission"
      self.job_id=self.Submit()
    else:
      self.job_id=None

  @property
  def sge_info(self):
      '''
      Info provided by the grid engine on our job.
      '''
      if self._sge_info:
          return self._sge_info
      self._sge_info = QacctInfo(self.job_id)
      return self._sge_info
      
  def _AssembleQSubCommand(self):
    '''
    Assembles the qsub command to submit this job. To avoid problems with 
    escaping, we put the arguments into a list that can then directly be passed
    to subprocess.Popen.
    '''
    print self.command
    print "start to assemble qsub command"

    qsub = [config.QSUB_BIN, '-b', 'y', '-S', '/bin/tcsh', '-N', self.name]
    if self.inherit_env:
      qsub.extend(['-V'])
    if self.cwd:
      qsub.append('-cwd')
    if self.hold_jid:
      qsub.append('-hold_jid')
      qsub.append(self.hold_jid)
    m_l_param = ['-l']
    print "look at the memory"
    if self.memory:
      # We just take memory for the WHOLE job but 'membycore' requires memory
      # per core. So we divide...
      if self.cpu and self.cpu > 1:
          membycore = float(self.memory[0]) / self.cpu
          if membycore < 0.001:
              membycore = 0.001
      else:
          membycore = self.memory[0]
      m_l_param.append('membycore=%.3f%s' % (membycore, self.memory[1]))
    print "look at stdout stuff"
    if len(m_l_param) > 1:
      qsub.extend(m_l_param)
    if self.stdout:
      if self.stdout==self.stderr or not self.stderr:
        qsub+=['-o', str(self.stdout), '-j', 'y']
      else:
        qsub+=['-o', str(self.stdout), '-e', str(self.stderr)]
    print "look at cpus"
    if self.cpu and self.cpu > 1:
      queue_name=self.queue.split('@@')[0]
      qsub+=['-pe', SGEJob.MT_ENVS[queue_name], str(self.cpu)]
    qsub+=['-q', self.queue]
    print "look whether its stringlike"
    print self.command
    if _IsStringLike(self.command):
      qsub+=shlex.split(self.command)
    else:
      qsub+=self.command

    print "returning the crap"
    return qsub

  def Submit(self):
    '''
    Submits the job to the queue and returns the job id
    '''
    qsub=self._AssembleQSubCommand()

    print "open subprocess"
    qsub_ps = subprocess.Popen(qsub, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    print "did it"
    pipe = Queue.Queue()
    # time is in seconds
    qsub_deadline = threading.Timer(300, _QsubTimeout, [qsub_ps, ' '.join(qsub),
                                                        pipe])
    qsub_deadline.start()
    j_re = re.compile(r'Your job (\d+) ')
    jobid = None
    for line in qsub_ps.stdout.readlines():
      m_jobid = re.search(j_re, line)
      if m_jobid:
        jobid = m_jobid.group(1)
        self.submitted=True
        break
    qsub_deadline.cancel()
    try:
       exc = pipe.get(block=False)
    except Queue.Empty:
      pass
    else:
      raise RuntimeError(exc)
    qsub_ps.wait()
    if qsub_ps.returncode != 0:
      err_out = ''
      for line in qsub_ps.stderr.readlines():
        err_out += line
      raise RuntimeError('qsub command "%s" did not exit with 0. Output:"\n%s"'\
                         % (' '.join(qsub), err_out))
    if not jobid:
      raise RuntimeError('Could not fetch job id for qsub command "%s"' % \
                         ' '.join(qsub))
    print "finished submission"

    return jobid


  def Wait(self):
    '''
    This method blocks until the job has finished
    '''
    while self.IsAlive():
      time.sleep(10)

  def IsAlive(self):
    '''
    Checks whether the job is still alive
    '''
    if not self.submitted:
      return False
    return IsAlive(self.job_id)

  def GotKilled(self):
      '''
      Check wether the grid engine killed our job.
      '''
      if not self.submitted:
          return False
      return self.sge_info.GotKilled()

  def GotKilledOfMem(self):
      '''
      Check wether the grid engine killed our job.
      '''
      if not self.submitted:
          return False
      return self.sge_info.GotKilledOfMem()

  def GotKilledOfTime(self):
      '''
      Check wether the grid engine killed our job.
      '''
      if not self.submitted:
          return False
      return self.sge_info.GotKilledOfTime()

  def Kill(self):
    '''
    Kills the job
    '''
    if not self.submitted:
      return
    subprocess.Popen([config.QDEL_BIN, self.job_id], shell=False)

class SGEJobList(list):
  '''
  Having a list of jobs with a convenience wait function.
  '''
  def Wait(self):
    '''
    Wait for the list of jobs to be finished.
    '''
    for job in self:
      job.Wait()

def WorkItems(filename, sh_lock, chunk_size=1):
  '''
  Simple function to facilitate distributing the workload to several processses. 
  The processes communicate via the filesystem, so it's even OK if the 
  processes run on different machines as long as they have access to the file 
  pointed to by filename.

  Each line in the file given by filename is ought to be one work item. As soon 
  as one process gets the work item, it is removed from the file. This means 
  that, if the process crashes during the execution, no other process will get 
  to see the work item.
  '''
  not_empty=True
  sh_lock_fd=open(sh_lock, 'w')
  while not_empty:
    fcntl.lockf(sh_lock_fd, fcntl.LOCK_EX)
    fd=open(filename, 'r+')
    lines=fd.readlines()
    fd.seek(0)
    fd.truncate()
    actual_chunk_size=min(len(lines), chunk_size)
    if actual_chunk_size==0:
      not_empty=False
      continue
    remaining=lines[actual_chunk_size:]
    fd.write(''.join(remaining))
    fd.flush()
    fcntl.lockf(sh_lock_fd, fcntl.LOCK_UN)
    for line in lines[:actual_chunk_size]:
      yield line.strip()

def ConcurrentWrite(filename, message, clear=False):
  '''
  Writes to file with given name. To avoid concurrent access of the file, 
  creates a lock file of the name filename.lock and make sure none else is 
  accessing it. If the file/the directory does not exist, it is created first.
  '''
  dirname=os.path.dirname(filename)
  if len(dirname) and not os.path.exists(dirname):
    os.makedirs(dirname)
  sh_lock='%s.lock' % filename
  sh_lock_fd=open(sh_lock, 'w')
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_EX)
  flags='a'
  if clear:
    flags='w'
  fd=open(filename, flags)
  fd.write('%s\n' % (message))
  fd.close()
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_UN)
  sh_lock_fd.close()

def ConcurrentRead(filename):
  '''
  Safely read a file which may be accessed (for reading/ writing) by other
  processes.

  :param filename: File to be read.
  :type filename: :class:`str`
  '''
  sh_lock='%s.lock' % filename
  if not os.path.exists(filename):
    return None
  sh_lock_fd=open(sh_lock, 'r')
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_SH)
  message=[l.strip() for l in open(filename, 'r')]
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_UN)
  sh_lock_fd.close()
  return message

def WriteStatus(filename, message, job_id=None, task_type=None):
  '''
  Write current job status to file
  '''
  if not job_id:
    job_id=ReadStatus(filename).split()[0].strip()
  elif task_type:
    task_type = task_type.upper()
    if task_type != 'SGE':
      task_type=socket.gethostname()
    job_id='%s@%s' % (str(job_id), task_type)
  msg='%s %s' % (job_id, message)
  ConcurrentWrite(filename, msg, clear=True)


def ReadStatus(filename):
  '''
  Read status from file.
  '''
  sh_lock='%s.lock' % filename
  if not os.path.exists(filename):
    return 'STATUS_FILE_NOT_FOUND_%s FAILED' % filename
  sh_lock_fd=open(sh_lock, 'r')
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_SH)
  message='\n'.join(open(filename, 'r').readlines())
  fcntl.lockf(sh_lock_fd, fcntl.LOCK_UN)
  sh_lock_fd.close()
  return message

def CheckStatus(status_file):
  '''
  Read status from file and does some additional sanity checks to make sure the 
  job hasn't failed.
  '''
  task_type='SGE'
  if not os.path.exists(status_file):
    return None
  
  line=ReadStatus(status_file).strip()
  job_id, status=line.split()
  tmp=job_id.split('@')
  if len(tmp) == 2:
    task_type=tmp[1]
    job_id=tmp[0]
  status=status.strip()
  if status in ('COMPLETED', 'FAILED'):
    return status
  if task_type == 'SGE':
    alive=IsAlive(job_id)
  else:
    if socket.gethostname() != task_type:
      return status
    alive=True
    try:
      os.kill(int(job_id), 0)
    except OSError:
      alive=False
  if alive:
    if status in ('QUEUEING', 'RUNNING'):
      return status
  else:
    if status in ('RUNNING', 'QUEUEING'):
      WriteStatus(status_file, 'FAILED', job_id)    
      return 'FAILED'
  raise RuntimeError('unknown job status %s' % status)

__all__=(
  'SGEJob', 'SGEJobList', 'ReadStatus', 'WriteStatus', 'ConcurrentWrite',
  'WorkItems', 'IsAlive', 'CheckStatus',
)
