# plasmo_genomics/utils/helpers.py

from pathlib import Path

def ensure_dir(path: Path):
    if not path.exists():
        print(f"📁 Output directory '{path}' does not exist. Creating it now...")
        path.mkdir(parents=True, exist_ok=True)
    else:
        print(f"📁 Output directory '{path}' already exists. Using it.")


import shutil
import sys

def check_tool_installed(tool_name):
    if shutil.which(tool_name) is None:
        print(f"❌ Required tool '{tool_name}' is not installed or not in your PATH.")
        sys.exit(1)
    else:
        print(f"🛠️ Found tool: {tool_name}")
