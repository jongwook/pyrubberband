from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

import sys
import imp
import glob

version = imp.load_source('pyrubberband.version', 'pyrubberband/version.py')

extra_compile_args = []
if 'win' in sys.platform:
    extra_compile_args += ['/DUSE_KISSFFT', '/DUSE_SPEEX']
else:
    extra_compile_args += ['-DUSE_KISSFFT', '-DUSE_SPEEX', '-DUSE_PTHREADS']

setup(
    name='pyrubberband',
    version=version.version,
    description='Python module to wrap rubberband',
    author='Brian McFee',
    author_email='brian.mcfee@nyu.edu',
    url='http://github.com/bmcfee/pyrubberband',
    download_url='http://github.com/bmcfee/pyrubberband/releases',
    packages=find_packages(),
    long_description="""A python module to wrap rubberband.""",
    classifiers=[
        "License :: OSI Approved :: ISC License (ISCL)",
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Multimedia :: Sound/Audio :: Analysis",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
    ],
    keywords='audio music sound',
    license='ISC',
    install_requires=[
        'six',
        'pysoundfile>=0.8.0',
        'Cython',
    ],
    extras_require={
        'docs': ['numpydoc'],
        'tests': [
            'pytest',
            'pytest-cov',
            'pytest-pep8'
        ]
    },
    test_require=[
        'pytest',
        'pytest-cov',
        'pytest-pep8'
    ],
    ext_modules=cythonize(Extension(
        'pyrubberband.rubberband',
        sources=['pyrubberband/rubberband.pyx'] + glob.glob('rubberband/**/*.c*'),
        include_dirs=['rubberband'],
        extra_compile_args=extra_compile_args,
        language='c++'
    ))
)
