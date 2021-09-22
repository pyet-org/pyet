from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# Get the version.
version = {}
with open("pyet/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='pyet',
    version=version['__version__'],
    url='https://github.com/phydrus/pyet',
    license='MIT License',
    author='Matevz Vremec, Raoul Collenteur',
    author_email='matevz.vremec@uni-graz.at, raoul.collenteur@uni-graz.at',
    description='pyet - Estimation of Potential Evaporation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    test_suite='tests',
    project_urls={
            'Source': 'https://github.com/phydrus/pyet',
            'Tracker': 'https://github.com/phydrus/pyet/issues',
            'Help': 'https://github.com/phydrus/pyet/discussions'
    },
    classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Other Audience',
            'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Topic :: Scientific/Engineering :: Hydrology',
    ],
    install_requires=['numpy>=1.15', 'matplotlib>=2.0', 'pandas>=1.0',
                      'scipy>=1.1'],
    packages=find_packages(exclude=[]),
)
