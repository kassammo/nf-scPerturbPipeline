# Author: Mohamed Kassam
# Date: 2026-07-19
from pathlib import Path
import subprocess, tempfile, csv

def test_help():
 r=subprocess.run(['python','bin/validate_inputs.py','--help'],capture_output=True,text=True)
 assert r.returncode==0
