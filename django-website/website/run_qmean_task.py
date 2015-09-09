import task
import os
import qmean_server_config as config

class RunQMEAN(task.Task):
  def __init__(self,project_dir):
    self.project_dir = project_dir
 
  def AssembleCommand(self):
    cmd = [config.OST_BIN, 
           os.path.join(config.QMEAN_SERVER_BASEDIR,"run_qmean.py"), 
           os.path.join(self.project_dir,"input"),
           os.path.join(self.project_dir,"output")]

  def CalculateMemConsumption(self):
    return (8.0,'G')

  task.Task.__init__(self, os.path.join(project_dir,"output"),
                           os.path.join(project_dir,"status"),
                           'run_qmean')   






