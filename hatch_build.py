from hatchling.builders.hooks.plugin.interface import BuildHookInterface
import sys, subprocess
class CustomBuildHook(BuildHookInterface):
    
    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r","requirements-backend.txt"])
