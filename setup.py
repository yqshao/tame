from setuptools import setup, find_packages

setup(name='tame',
      version='dev',
      description='Trajectory Analysis Made Easy',
      url='https://github.com/yqshao/mdppp',
      author='Yunqi Shao',
      author_email='yunqi_shao@yahoo.com',
      license='BSD',
      packages=find_packages(),
      install_requires=['numpy'],
      entry_points={
          'console_scripts': ['tame=tame.recipes.bin:main']
      }
)
