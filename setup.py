import sys
import setuptools

if sys.version_info[0] < 3:
    sys.stdout.write('package requires python3')
    sys.exit(1)

setuptools.setup()
