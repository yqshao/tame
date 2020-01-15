from setuptools import setup, find_packages

setup(name='mdppp',
      version='dev',
      description='Molecular Dynamics Post-Processing Programm',
      url='https://github.com/yqshao/mdppp',
      author='Yunqi Shao',
      author_email='yunqi_shao@yahoo.com',
      license='BSD',
      packages=find_packages(),
      install_requires=['numpy'],
      entry_points={
          'console_scripts': ['mdppp=mdppp.recipes.bin:main']
      }
)
