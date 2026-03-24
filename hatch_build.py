from hatchling.builders.hooks.plugin.interface import BuildHookInterface
import sys, subprocess
# import logging
# logger = logging.getLogger(__name__)
class CustomBuildHook(BuildHookInterface):
    def initialize(self,*param,**kwargs):
        with open("./analysator/miscellaneous/_commithash.py","w") as file:
            commithash=subprocess.run(["git","rev-parse","HEAD"],capture_output=True).stdout.decode('utf-8')
            file.write(f"commithash='{commithash.strip('\n')}'\n")
            file.close()
       
    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str)->None:
        out=subprocess.run([sys.executable, "-m", "pip", "install", "-r","requirements-backend.txt","-v"],capture_output=True,check=True)
        #logger.error("TESTINGTESTING ALERTA"+out.stdout.decode('utf-8')) #Doesnt work
        #print(out) #Also doesnt work
        # self.app.display_error(out.stdout.decode('utf-8')) #Does not work
