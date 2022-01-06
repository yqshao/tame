import re
from setuptools import setup, find_packages

with open('tame/__init__.py') as f:
    version = re.search("__version__ = '(.*)'", f.read()).group(1)

setup(name='tame',
      version=version,
      description='Trajectory Analysis Made Easy',
      url='https://github.com/yqshao/mdppp',
      author='Yunqi Shao',
      author_email='yunqi_shao@yahoo.com',
      license='BSD',
      packages=find_packages(),
      install_requires=['numpy>=1.8', 'click>=7.0'],
      entry_points={
          'console_scripts': ['tame=tame.recipes.bin:main']
      }
)
