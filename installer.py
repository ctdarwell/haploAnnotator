import importlib, subprocess, sys

packages = sys.argv[1:]
print('Installing', packages)
[subprocess.check_call(['pip', 'install', pkg]) for pkg in packages if not importlib.util.find_spec(pkg)]



