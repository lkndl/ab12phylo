import sys

from setuptools import setup

if sys.version_info[0] < 3:
    sys.stdout.write('package requires python3')
    sys.exit(1)

from ab12phylo.__init__ import __version__, __author__, __email__, __license__

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
      packages=['ab12phylo', 'ab12phylo_gui'],
      package_dir={'ab12phylo': 'ab12phylo', 'ab12phylo_gui': 'ab12phylo_gui'},
      include_package_data=True,
      zip_safe=True,
      entry_points={'console_scripts': ['ab12phylo = ab12phylo.main:_main',
                                        'ab12phylo-visualize = ab12phylo.phylo:_visualize',
                                        'ab12phylo-view = ab12phylo.phylo:_view',
                                        'ab12phylo-add-xml = ab12phylo.blast:_add_xml'],
                    'gui_scripts': ['ab12phylo-gui = ab12phylo_gui.gtk_main:_main']},
      install_requires=['biopython', 'pyyaml', 'jinja2', 'lxml', 'xmltramp2', 'toyplot',
                        'toytree<=1.2.0', 'numpy', 'pandas', 'matplotlib', 'svgutils'],
      classifiers=['Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent'],
      keywords=['bioinformatics', 'phylogenetics', 'population genetics'],
      python_requires='>=3.6')
