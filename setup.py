#!/usr/bin/env python
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    # Metadata
    name='WAVEPAL',
    version='3.3',
    author='Guillaume Lenoir',
    author_email='guillaume.lenoir@hotmail.com',
    url='http://www.climate.be/u/lenoir',
    description='WAVEPAL Python project',
    long_description=long_description,
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: MacOS X',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    keywords='WAVEPAL setuptools development',
    packages=find_packages(exclude=['carmcmc', 'test', 'examples']),
    install_requires=['numpy', 'scipy', 'tqdm', 'carmcmc', 'matplotlib', 'acor'],
    zip_safe=True,
)
