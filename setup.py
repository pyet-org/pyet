from setuptools import setup

setup(
    name='pyet',
    version='0.2b',
    packages=[''],
    url='https://github.com/phydrus/pyet',
    license='GNU General public version 3.0',
    author='Matevz Vremec, Raoul Collenteur',
    author_email='matevz.vremec@uni-graz.at, raoul.collenteur@uni-graz.at',
    description='pyet',
    test_suite='tests',
    install_requires=['numpy>=1.15', 'matplotlib>=3.0', 'pandas>=1.0',
                      'scipy>=1.1'],
)
