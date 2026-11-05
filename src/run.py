import sys
import subprocess
from pathlib import Path

def run():
    subprocess.call([sys.executable, f"{Path(__file__).parent}/AIMS.py", *sys.argv[1:]])

if __name__ == "__main__":
    run()
