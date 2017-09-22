Classification
==============

Kpop population objects have builtin tools for building classifiers from
population objects. As with any classification task, you must provide a labeled
training data set and the classifier algorithm will be trained to replicate
those labels and generalize to new sample points. A very basic classification
task can start with a population object and a list of labels:

>>> from kpop import Population
>>> pop = Population.random(5, 10)
>>> cls = pop.classification(['A', 'A', 'B', 'B', 'A'])

The method returns a trained classifier object that associates each individual
in the population with the given labels. Notice that we created a random
population with 5 individuals and we had to provide the same number of labels.


The classifier implements
several different classification algorithms that can be accessed either
using the ``Population.classification(<method>)`` method or using the
corresponding attribute ``Population.classification.<method>()``.




Api docs
--------

.. autoclass:: kpop.population.classification
    :members:
