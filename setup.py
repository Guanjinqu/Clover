from setuptools import setup

setup(name='dna-clover',
      version='1.2',
      description='Tree structure-based efficient DNA clustering for DNA-based data storage',
      url='http://github.com/storborg/funniest',
      author='Guanjin qu',
      author_email='guanjinqu@tju.edu.cn',
      license='GNU',
      packages=['clover','tests'],
      zip_safe=False,
      entry_points = {
        'console_scripts': ['clover = clover.main:main'],
    })