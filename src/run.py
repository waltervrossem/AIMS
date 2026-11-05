import sys
import shutil
import subprocess
from pathlib import Path

this_dir = Path(__file__).parent

def run():
    subprocess.call([sys.executable, f"{this_dir}/AIMS.py", *sys.argv[1:]])


def get_default_config():
    shutil.copy2(f'{this_dir}/AIMS_configure.py', f'./AIMS_configure.py')
