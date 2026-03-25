import sys, subprocess
from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from packaging.version import Version

class CustomBuildHook(BuildHookInterface):
    def initialize(self,*param,**kwargs):
        with open("./analysator/miscellaneous/_commithash.py","w") as file:
            commithash=subprocess.run(["git","rev-parse","HEAD"],capture_output=True).stdout.decode('utf-8')
            commithash=commithash.strip("\n")
            file.write(f"commithash='{commithash}'"+"\n")
            file.close()

