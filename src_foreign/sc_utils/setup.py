from setuptools import setup

setup(name='sc_utils',
      version='0.1',
      description='Single cell analysis',
      url='http://github.com',
      author='Kai Zhang',
      author_email='kai@kzhang.org',
      license='MIT',
      packages=['sc_utils'],
      entry_points = {
        'console_scripts': ['sc_utils=sc_utils.command_line:main'],
      },
      install_requires=[
          'gensim',
          'scipy',
          'numpy',
          'python-igraph',
          'umap-learn',
          'leidenalg'
      ],
      zip_safe=False)