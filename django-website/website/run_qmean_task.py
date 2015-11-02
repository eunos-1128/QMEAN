import task
import os
from sm import config
from django.conf import settings

class RunQMEAN(task.Task):
  def __init__(self,project_dir):
    self.project_dir = project_dir
 
    task.Task.__init__(self, os.path.join(project_dir,"output"),
                             os.path.join(project_dir,"status"),
                             'run_qmean')

  def AssembleCommand(self):
    cmd = [os.path.join(config.SM_ROOT_DIR,"bin","sm"), 
           os.path.join(settings.BASE_DIR,"website","run_qmean.py"), 
           os.path.join(self.project_dir,"input"),
           os.path.join(self.project_dir,"output"),
           settings.CACHE_DIR]
           
    return cmd

  def CalculateMemConsumption(self):
    return (8.0,'G')

   






