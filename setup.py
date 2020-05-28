import sys
from setuptools import setup

if sys.version_info[0] < 3:
    sys.stdout.write('package requires python3')
    sys.exit(1)

setup(name='ab12phylo',
      version='0.1a.2',
      author='Leo Kaindl',
      author_email='leo.kaindl@tum.de',
      license='MIT',
      description='An integrated pipeline for maximum likelihood '
                  'phylogenetic tree inference from ABI sequencing data',
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      url='https://gitlab.lrz.de/leokaindl/ab12phylo',
      packages=['ab12phylo'],
      package_dir={'ab12phylo': 'ab12phylo'},
      include_package_data=True,
      zip_safe=True,
      data_files=[('docs', ['docs/RAXML.md'])],
      entry_points={'console_scripts': ['ab12phylo = ab12phylo.main:_main',
                                        'ab12phylo-visualize = ab12phylo.phylo:_visualize',
                                        'ab12phylo-view = ab12phylo.phylo:_view']},
      install_requires=['biopython', 'dendropy', 'pyyaml', 'jinja2',
                        'toyplot', 'toytree', 'numpy', 'pandas'],
      classifiers=['Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent'],
      keywords=['bioinformatics', 'phylogenetics'],
      python_requires='>=3.8')
