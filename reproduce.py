#!/usr/bin/env python3
"""Rebuild the certificate and verify it from scratch, offline, with Python 3."""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

root = Path(__file__).resolve().parent
start = time.perf_counter()
subprocess.run([sys.executable, str(root/"build_certificate.py")], check=True)
subprocess.run([sys.executable, str(root/"verify.py"), "--output", str(root/"verification.json")], check=True)
subprocess.run([sys.executable, str(root/"audit_checker.py")], check=True)
summary = {"reproduction_passed": True,
           "wall_seconds": time.perf_counter()-start,
           "python": sys.version.split()[0],
           "files_sha256": {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                            for name in ("build_certificate.py", "verify.py", "audit_checker.py", "certificate.json")}}
(root/"reproduction.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps(summary, indent=2))
