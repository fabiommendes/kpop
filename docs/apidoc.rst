=============
API Reference
=============

API documentation for the :mod:`kpop` module.


Individual
==========

Each element of a population is an instance of :class:`kpop.Individual`. An
:class:`kpop.Individual` behave similarly as a list of genetypes or as an
2D array of genotypes.

.. autoclass:: kpop.Individual
   :members:
   :inherited-members:


Population objects
==================

The main type in the kpop package is :class:`kpop.Population`. A population is basically
a list of individuals. It has a similar interface as a Python's list or a Numpy
array.

.. autoclass:: kpop.Population
   :members:
   :inherited-members:


Population vs Multipopulation
-----------------------------

Kpop uses two classes to represent populations that have basically the same
interface. A MultiPopulation is basically a population structured with many
sub-populations.

.. autoclass:: kpop.MultiPopulation
   :members:
   :inherited-members:


The ``.plot`` attribute
-----------------------

Each :class:`kpop.Population` or :class:`kpop.MultiPopulation` instance have a
``.plot`` attribute that defines a namespace with many different plotting
utilities.

.. autoclass:: kpop.population.attr_plot.PlotAttribute
   :members:


Other utility types
===================

Representing probabilities
--------------------------

.. autoclass:: kpop.prob.Prob
   :members:


Utility modules
===============

Plotting
--------

:mod:`kpop.plots` contains a few useful plotting functions based on
matplotlib.

.. automodule:: kpop.plots
   :members:


Loading objects
---------------

Functions from the :mod:`kpop.loaders` module are responsible for loading
Population objects from files.

.. automodule:: kpop.loaders
   :members:



