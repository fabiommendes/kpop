.. image:: https://readthedocs.org/projects/kpop/badge/?version=latest
   :target: https://kpop.readthedocs.io/en/latest/
   :alt: Read The Docs

.. image:: https://codecov.io/gh/fabiommendes/kpop/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/fabiommendes/kpop
   :alt: Codecov

.. image:: https://travis-ci.org/fabiommendes/kpop.svg?branch=master
   :target: https://travis-ci.org/fabiommendes/kpop
   :alt: Travis-CI

.. image:: https://codeclimate.com/github/fabiommendes/kpop/badges/gpa.svg
   :target: https://codeclimate.com/github/fabiommendes/kpop
   :alt: Code Climate

.. image:: https://www.quantifiedcode.com/api/v1/project/989e3bf4fb4f4afe9fdf762d819b3fe4/badge.svg
  :target: https://www.quantifiedcode.com/app/project/989e3bf4fb4f4afe9fdf762d819b3fe4
  :alt: Code issues


Kpop is a Python package to perform population genetics analysis that
integrates traditional methods such as PCA and LDA (Latent Dirichilet Analysis,
a.k.a the algorithm introduced by Pritchard in the Structure program) with
machine learning.

Command line interface
----------------------

Kpop can be used either as a library or as a command line interface application.
Although the former is more flexible and powerful, some simple transformations
and analysis can be done directly on the command line without touching any Python
code.

Conversion from kpop to other common formats. Export to Structure database and
param files::

    $ kpop export pop.csv -f structure


Simple graphics. Show the result of a PCA analysis::

    $ kpop show pop.csv -m pca

Clusterization/sharp estimate of parental affliations::

    $ kpop cluster pop.csv -m kmeans -k 2

Soft cluster/computation of admixture coefficients::

    $ kpop admix pop.csv -m kmeans -k 2

Basic statistics::

    $ kpop stats pop.csv

Creation of random synthetic populations::

    $ kpop create pop.csv --size 100 --num-loci 200

Start shell loading all of kpop symbols and all population files and in the
current directory::

    $ kpop shell

We encourage you to explore more options by using the builtin help utility::

    $ kpop --help           # global help
    $ kpop show --help      # help for a specific command


Kpop data formats
-----------------

Kpop defines two file formats that represents populations. Pickle is a binary
format that can be used to store the full state of a population, including any
additional fields and meta data you may have created. Pickle is used internally
by Python to serialize objects and is the fastest and most flexible format.

A major drawback of using Pickle is that other programs will not understand it.
It can also change across Python versions, so if you create a database with a
later Python version, it might not work when loaded from older Python
interpreters. Pickle may also break after a major version upgrade of Kpop
itself.

Pickle is fast and convenient, but it is not an archival and data exchange
format. For that, Kpop uses simple CSV files. CSV is limited and is not the
most efficient format both in terms of loading speed and disk usage. It is
however easy to produce and you can even use your favorite spreadsheet
program to create/edit an CSV file.

Kpop expects that CSV files should have a certain structure. The :func:`kpop.load_csv`
function can adapt to different formats, but the command line interface expects
a more or less rigid configuration. Your CSV file must have a single header line
for which Kpop understands a few column names:

label:
    label for a single individual.
index:
    a numeric index. Kpop will ignore this column if label is given or otherwise
    use it as a label.
population:
    a label or index for the population that each individual belongs to.
gender:
    arbitrary gender label string. Not restricted to male/female.
age:
    a numeric (float) value representing age. You decide if this number means
    years, days, minutes, simulation ticks, etc.
phenotype:
    arbitrary string describing phenotype.
meta information:
    any column named as "#some-name" will be treated as arbitrary meta
    information attached to each individual. This data is stored, but does not
    influence any analysis performed by Kpop.

All other columns are treated as genetic marker names. The content of each marker
is a string of N characters in which each character represents an specific
allele. Kpop prefers numeric identifiers (e.g.: 12, 11, 22) vs letters (eg.: aA,
AB, aa, etc), but it also accepts the later. By default, it treats ``0``, the
dash character (``-``) and empty cells as missing data.

A typical Kpop CSV file will be like the following::

    label,population,rs123,rs1234,rs42
    john,uk,12,02,11
    paul,uk,22,22,
    psy,korea,--,21,22

In example above, "john" has a missing allele in the second locus and "paul" and
"psy" have no data for an entire locus (the third and the first, respectively).


Kpop python interface
---------------------

In order to enjoy the full power of Kpop, it is necessary to use it from Python.
Kpop can be used as a library by importing it in python code:

.. code-block:: python

    import kpop

    pop = kpop.Population.random(10, 100)
    ...

If you are just exploring, it might be more useful to just open the Python shell
or a Jupyter notebook using one of the commands::

    $ kpop shell
    $ kpop shell --notebook

It will start a Jupyter shell (or notebook) that already loads all symbols in
the Kpop namespace and

Users are refered to the :doc:`API Reference<apidoc>`

