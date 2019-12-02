from setuptools import setup
setup(
    name='etmodule',
    version='0.1',
    packages=[''],
    url='https://github.com/matevzvremec/etmodule',
    license='GNU General public version 3.0',
    author='Matevz Vremec',
    author_email='matevz.vremec@uni-graz.at',
    description='etmodule'
    install_requires=['numpy>=1.15', 'matplotlib>=2.0', 'pandas>=0.23',
                      'scipy>=1.1']
)
