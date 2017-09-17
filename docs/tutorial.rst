========
Tutorial
========

Getting started
===============

We are assuming you already installed Kpop. The easiest way to start is to
simply type the following command on the terminal::

    $ kpop shell 

This will start a Python session with all basic Kpop functionality available. 
It also also tries to load all population files in the current directory, so
they become easily available.

The kpop shell is just a convenience method of starting a IPython shell with a
few useful imports:

.. code-block: python

    from kpop import *
    import numpy as np
    import scipy as sp
    import matplotlib.pyplot as plt

This is useful for interactive and exploratory work. However, for more serious
jobs you probably should make those imports manually and avoid the start import
in the first line.


Basic Kpop concepts
===================

You notice that most of Kpop interactions go through two main object types:
:class:`kpop.Individual` and :class:`kpop.Population`. Let us start with the first
of these two, :class:`kpop.Individual`, which represent single individuals by
their corresponding genotypes. 

Individual
----------

An :class:`kpop.Individual` instance behaves basically as a list of genotype
values. Kpop represents genotypes by numbers, where zero is used to encode missing
data and numbers above one represent each allele. We can start a new 
individual by constructing it from a list of pairs of numbers:

>>> ind = Individual([[1, 1], [1, 2], [2, 2], [1, 2]])

This is a genotype with 4 loci of biallelic data. You might expect it behave
just as a list of genotypes for each locus. It accepts Python indexing, slicing
and iteration:

>>> ind[0]
array([1, 1], dtype=uint8)

>>> [(1 in locus) for locus in ind]
[True, True, False, True]

:class:`kpop.Individual` objects can also be inspected in several ways.

>>> ind.num_loci, ind.ploidy, ind.is_biallelic
(4, 2, True)

You should use the autocomplete feature of Kpop's shell to discover more 
attributes. Just type ``ind.`` and hit the <tab> key to see a list of 
completions. Some of those options are methods (you will notice it by the 
open-close parens at the end of their names). In order to get help on the 
methods behavior and signature, just use the ``?`` helper as bellow

>>> ind.breed?                                                  # doctest: +SKIP
Signature: ind.breed(other, id=None, **kwargs)
   doctest:
Breeds with other individual.
<NEWLINE>
Creates a new genotype in which features are selected from both
parents.
File:      ~/git/bio/kpop/src/kpop/individual.py
Type:      method

You will notice that if you print an Individual in the terminal it will shown
with the following notation

>>> ind
Individual('ind: 11 12 22 12')
 
This is actually a different way to construct :class:`kpop.Individual` instances.
The first part in the string before the column is a label used to identify the
given individual and everything on the right hand side is its genotype. 

Let us create a second individual to interact with the first.

>>> ind2 = Individual('ind2: 22 11 12 12')
>>> ind2.breed(ind)                                             # doctest: +SKIP
Individual('ind2_: 21 12 12 12')

Of course, handling a handful of individuals is not very useful. Let us create a
list of individuals by drawing samples from an specific probability. First, 
define a list of probabilities for each allele in each loci

>>> freqs = [[0.1, 0.9], [0.5, 0.5], [0.9, 0.1], [0.5, 0.5]]

Now we can create a random individual using the ``from_freqs`` method of the 
Individual class

>>> random_ind = Individual.from_freqs(freqs)

... and now we create a bunch:

>>> list_of_individuals = []
>>> for _ in range(10):
...     new_ind = Individual.from_freqs(freqs)
...     list_of_individuals.append(new_ind) 


Population
----------

Now that we have a bunch of individuals, we can make a population. Of course
we could use the list of individuals directly, but Kpop provides the much more
convenient :class:`kpop.Population` type to represent a group of individuals. 

>>> popA = Population(list_of_individuals, id='A')
>>> popA                                                        # doctest: +SKIP
  ind1: 22 21 12 22
  ind2: 22 11 11 21
  ind3: 22 11 11 21
  ind4: 22 11 11 21
  ind5: 22 11 11 12
  ind6: 22 22 11 21
  ind7: 22 11 11 21
  ind8: 22 21 11 22
  ind9: 22 12 11 12
 ind10: 22 22 11 21

We created the Population object from a list of individuals and gave it an 
optional label. The label is used to identify the population in several different
contexts such as clustering, plotting, etc.

Just like :class:`kpop.Individual` instances, :class:`kpop.Population` objects
have many associated methods and attributes. You can explore it by typing
``popA.`` and hitting the <tab> key (you will notice it is way more complex than
Individual instances).

In population genetics we are usually interested in comparing different
populations rather than different individuals in the same population. We can
easily create a new random population using the Population.make_random
function:

>>> popB = Population.random(10, num_loci=4, id='B')

This will create a new population with 10 individuals and 4 loci. Now, let us
compose this population with the previous one by creating a new generation that
breeds individuals from the first population with the second

>>> popC = popA.simulation.breed(popB, size=15, id='C')

We can combine all sub-populations into a single population containing all
individuals by simply adding the population objects together

>>> pop_all = popA + popB + popC


This creates a :class:`kpop.MultiPopulation` object which behaves essentially as
a Population, but keeps track of sub-structuring. 


Visualization
=============

Kpop implements a few visualization methods through the Population.plot 
attribute. The ``population.plot.?`` namespace has methods for dimensionality
reduction (such as PCA), 


Statistics
==========





Admixture
=========

Admixture analysis is the task of estimating the admixture coefficients of each
individual in a population. This is the main concern of programs such as 
Structure and ADMIXTURE.


Projections
===========

All dimensionality reduction methods from the above section are implemented in
the ``population.projection`` namespace. Those methods provide the raw data for
dimensionality reduction and may be useful in contexts other than data 
visualization.

# TODO.


Clusterization
==========

Clusterization is the task of spliting data into separate groups without providing
a training set on correct classifications. This is often refered as "unsupervised
learning". Notice here that "unsupervised" does not mean "completely independent
of human intervention" since almost all clustering algorithms requires some
sort of tuning.

Kpop provides a few methods for performing clustering of individuals. They are
all implemented under the ``population.cluster`` namespace. 

# TODO.


Classification
==============

Differently from clustering, a supervised classification task learns from a 
dataset in which all items are classified with a corresponding label. A 
classification task is useful when it can generalize this mapping to data points 
outside of the training set.

In Population genetics this often maps to the sittuation in which we have a 
group of individuals with known parental populations and we want to classify 
additional specimens into one of those populations. Notice it is different from 
admixture analysis that tries to infer the fractions of DNA belonging to each 
parental population. Here the classification is sharp: the individual is said
to belong to a single parental population.

All classification methods live under the ``population.classification`` 
namespace.

# TODO

