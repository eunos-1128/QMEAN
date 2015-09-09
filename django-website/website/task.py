import sge
import os, time
import sys
import subprocess

class Task:
  """
  A task to be run on the SGE queueing system
  
  """

  def __init__(self, status_dir, status_file, task_name, wait=False, submit=True, use='SGE',
               cpu=None, queue='long.q'):

    self.task_type=use.upper()
    self.cpu = cpu
    self.status_file=status_file
    self.command=self.AssembleCommand()
    stdout = os.path.join(status_dir,
                          '%s.stdout' % task_name)
    stderr = os.path.join(status_dir,
                          '%s.stderr' % task_name)
    job_id=-1
    if self.task_type == 'SGE':
      try:
        print "start to create sge job"
        task_mem = self.CalculateMemConsumption()
        self.job=sge.SGEJob(self.command, submit=submit, inherit_env=True,
                            stdout=stdout, stderr=stderr, cpu=cpu, queue=queue,
                            name=task_name, memory=task_mem)
        print "finished sge job creation"
        if submit:
          job_id = self.job.job_id

      except Exception, e:
        sge.WriteStatus(self.status_file, 'FAILED', 0, task_type=use)
        raise e
          
    elif self.task_type == 'PROCESS':
      if submit:
        if not hasattr(stdout, 'write'):
          stdout = open(stdout, 'w')
        if not hasattr(stderr, 'write'):
          stderr = open(stderr, 'w')
        self.job = subprocess.Popen(self.command, stdout=stdout, stderr=stderr)
        job_id=self.job.pid

    if submit:
      sge.WriteStatus(self.status_file, 'QUEUEING', job_id, task_type=use)

    if submit and wait:
      if self.task_type == 'SGE':
        self.job.Wait()
      elif self.task_type == 'PROCESS':
        self.job.wait()
        sge.CheckStatus(self.status_file)
  
  def AssembleCommand(self):
    """
    Subclasses are required to implement this method and return the command
    line arguments to execute the task as a list of strings. This list will
    then be passed to either subprocess.Popen, or when using the SGE queuing
    system, the :class:`~sm.core.sge.SGEJob` class.

    :raises: :exc:`NotImplemented`
    """
    raise NotImplemented('the method needs to be implemented in subclasses')

  def CalculateMemConsumption(self):
    """
    Subclasses are required to implement this method and return the supposed
    memory consumption as a tuple like (4.0, 'G'). This string will then be
    passed to the SGE queueing system, the :class:`~sm.core.sge.SGEJob` class,
    if the task is called with use='SGE'.

    :raises: :exc:`NotImplemented`
    """
    raise NotImplemented('the method needs to be implemented in subclasses')

  @property
  def completed(self):
    """
    Whether the task has completed
    """
    return self.status=='COMPLETED'
  @property
  def failed(self):
    """
    Whether the task has failed
    """
    return self.status=='FAILED'
  @property
  def running(self):
    """
    Whether the task is still running
    """
    return self.status=='RUNNING'
  
  @property 
  def queueing(self):
    """
    Whether the task is queueing
    """

    return self.status=='QUEUEING'

  @property
  def is_alive(self):
    """
    Whether the task is alive, e.g. it is running or queueing.
    """

    return self.status in ('QUEUEING', 'RUNNING')

  @property
  def status(self):
    if self.task_type == 'PROCESS':
      self.job.poll()
    return sge.CheckStatus(self.status_file)


class TaskPool:
  """
  A task pool simplifies the process of running several tasks in 
  parallel, but keeps the number at bay. New tasks may be added 
  with :meth:`Add`. When the maximum number of tasks is reached,
  """
  def __init__(self, max_tasks):
    self.max_tasks = max_tasks
    self.tasks = []
    self.failed = []
    self.completed = []
  def WaitOne(self):
    while True:
      for i, task in enumerate(self.tasks):
        if not task.is_alive:
          del self.tasks[i]
          if task.failed:
            self.failed.append(task)
          else:
            self.completed.append(task)
          return
      time.sleep(1)

  def Wait(self):
    while len(self.tasks)>0:
      self.WaitOne()
  def Add(self, task, message, *args, **kwargs):
    while len(self.tasks) >= self.max_tasks:
      self.WaitOne()
    print message
    self.tasks.append(task(*args, **kwargs))

