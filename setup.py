from setuptools import setup, find_packages

setup(name='mdppp',
      version='dev',
      description='Molecular Dynamics Post-Processing Programm',
      url='https://github.com/yqshao/mdppp',
      author='Yunqi Shao',
      author_email='yunqi_shao@yahoo.com',
      license='BSD',
      packages=find_packages(),
      entry_points={
          'console_scripts': ['mdppp=mdppp.recipe.parser:main']
      }
)
