import sys
from setuptools import setup

if sys.version_info[0] < 3:
    sys.stdout.write('package requires python3')
    sys.exit(1)

__author__ = 'Leo Kaindl'
__email__ = 'leo.kaindl@tum.de'
__version__ = '0.2b.6'
__date__ = '19 August 2020'
__license__ = 'MIT'
__status__ = 'Beta'

setup(name='ab12phylo',
      version=__version__,
      author=__author__,
      author_email=__email__,
      license=__license__,
      description='An integrated pipeline for maximum likelihood '
                  'phylogenetic tree inference from ABI sequencing data',
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      url='https://gitlab.lrz.de/leokaindl/ab12phylo',
      packages=['ab12phylo'],
      package_dir={'ab12phylo': 'ab12phylo'},
      include_package_data=True,
      zip_safe=True,
      entry_points={'console_scripts': ['ab12phylo = ab12phylo.main:_main',
                                        'ab12phylo-visualize = ab12phylo.phylo:_visualize',
                                        'ab12phylo-view = ab12phylo.phylo:_view',
                                        'ab12phylo-add-xml = ab12phylo.blast:_add_xml']},
      install_requires=['biopython', 'pyyaml', 'jinja2', 'lxml', 'xmltramp2',
                        'toyplot', 'toytree', 'numpy', 'pandas'],
      classifiers=['Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent'],
      keywords=['bioinformatics', 'phylogenetics'],
      python_requires='>=3.8')
