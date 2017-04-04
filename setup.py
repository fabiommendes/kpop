import os
import codecs
from setuptools import setup, find_packages
from distutils.extension import Extension

try:
    if 'BUILDC' in os.environ:
        raise ImportError
    from Cython.Distutils import build_ext
    mod_ext = '.pyx'
except ImportError:
    mod_ext = '.c'
    build_ext = None

# Save version and author to __meta__.py
version = open('VERSION').read().strip()
dirname = os.path.dirname(__file__)
path = os.path.join(dirname, 'src', 'kpop', '__meta__.py')
meta = '''# Automatically created. Please do not edit.
__version__ = '%s'
__author__ = 'F\\xe1bio Mac\\xeado Mendes'
''' % version
with open(path, 'w') as F:
    F.write(meta)

setup(
    # Basic info
    name='kpop',
    version=version,
    author='Fábio Macêdo Mendes',
    author_email='fabiomacedomendes@gmail.com',
    url='https://github.com/fabiommendes/kpop/',
    description='A toolset for population genetics with bindings to common '
                'programs and algorithms.',
    long_description=codecs.open('README.rst', 'rb', 'utf8').read(),

    # Classifiers (see https://pypi.python.org/pypi?%3Aaction=list_classifiers)
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries',
    ],

    # Entry points
    entry_points={
        'console_scripts': [
            'kpop = kpop.__main__:main'
        ],
    },

    # Packages and dependencies
    package_dir={'': 'src'},
    packages=find_packages('src'),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'sklearn',
        'maputil',
        'sequtil',
        'pandas',
        'pyplink',
        'click',
        'colorama',
    ],
    extras_require={
        'dev': [
            'python-boilerplate[dev]',
        ],
    },
    cmdclass={"build_ext": build_ext} if build_ext else {},
    ext_modules=[
        Extension("kpop.admixture.linalg",
                  ["src/kpop/admixture/linalg" + mod_ext],
                  libraries=["m"],
                  include_dirs=['src'],
        ),
        Extension("kpop.admixture.util",
                  ["src/kpop/admixture/util" + mod_ext],
                  libraries=["m"],
                  include_dirs=['src'],
        ),
        Extension("kpop.admixture.objective",
                  ["src/kpop/admixture/objective" + mod_ext],
                  libraries=["m"],
                  include_dirs=['src'],
        ),
        Extension("kpop.admixture.likelihood",
                  ["src/kpop/admixture/likelihood" + mod_ext],
                  libraries=["m"],
                  include_dirs=['src'],
        ),
        Extension("kpop.admixture.em",
                  ["src/kpop/admixture/em" + mod_ext],
                  libraries=["m"],
                  include_dirs=['src'],
        ),
    ],

    # Other configurations
    zip_safe=False,
    platforms='any',
)