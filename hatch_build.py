from hatchling.builders.hooks.plugin.interface import BuildHookInterface
import sys, subprocess
from packaging.version import Version
# import logging
# logger = logging.getLogger(__name__)
import platform
if Version(platform.python_version()) <= Version("3.7"):
    quit()
class CustomBuildHook(BuildHookInterface):
    def initialize(self,*param,**kwargs):
        with open("./analysator/miscellaneous/_commithash.py","w") as file:
            commithash=subprocess.run(["git","rev-parse","HEAD"],capture_output=True).stdout.decode('utf-8')
            commithash=commithash.strip("\n")
            file.write(f"commithash='{commithash}'"+"\n")
            file.close()
       
    def finalize(self, *param,**kwargs)->None:
        out=subprocess.run([sys.executable, "-m", "pip", "install", "-r","requirements-backend.txt","-v"],capture_output=True,check=True)
        #logger.error("TESTINGTESTING ALERTA"+out.stdout.decode('utf-8')) #Doesnt work
        #print(out) #Also doesnt work
        # self.app.display_error(out.stdout.decode('utf-8')) #Does not work
