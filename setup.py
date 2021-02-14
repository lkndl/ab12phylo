import sys
from types import SimpleNamespace

from setuptools import setup
from setuptools.config import read_configuration

if sys.version_info[0] < 3:
    sys.stdout.write('package requires python3')
    sys.exit(1)

pkg_meta = SimpleNamespace(**read_configuration('setup.cfg')['metadata'])

setup(name='ab12phylo-lkndl',
      version=pkg_meta.version,
      author=pkg_meta.author,
      author_email=pkg_meta.author_email,
      license=pkg_meta.license,
      description='An integrated pipeline for maximum likelihood '
                  'phylogenetic tree inference from ABI sequencing data',
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      url='https://github.com/lkndl/ab12phylo',
      packages=['ab12phylo_cmd', 'ab12phylo'],
      package_dir={'ab12phylo_cmd': 'ab12phylo_cmd', 'ab12phylo': 'ab12phylo'},
      include_package_data=True,
      zip_safe=True,
      entry_points={'console_scripts': ['ab12phylo-cmd = ab12phylo_cmd.main:_main',
                                        'ab12phylo-visualize = ab12phylo_cmd.phylo:_visualize',
                                        'ab12phylo-view = ab12phylo_cmd.phylo:_view',
                                        'ab12phylo-add-xml = ab12phylo_cmd.blast:_add_xml'],
                    'gui_scripts': ['ab12phylo = ab12phylo.gtk_app:main']},
      install_requires=['biopython', 'pyyaml', 'jinja2', 'lxml', 'xmltramp2',
                        'toyplot', 'toytree<=1.2.0', 'numpy', 'pandas',
                        'matplotlib', 'svgutils', 'pillow', 'requests'],
      classifiers=['Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Operating System :: OS Independent'],
      keywords=['bioinformatics', 'phylogenetics', 'population genetics'],
      python_requires='>=3.6')
