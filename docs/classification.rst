Classification
==============

Kpop population objects have builtin tools for building classifiers from
population objects. As with any classification task, you must provide a labeled
training data set and the classifier algorithm will be trained to replicate
those labels and generalize to new sample points. A very basic classification
task can start with a population object and a list of labels:

>>> from kpop import Population
>>> pop = Population.random(5, 10)
>>> classifier = pop.classification(['A', 'A', 'B', 'B', 'A'])

The method returns a trained classifier object that associates each individual
in the population with the given labels. Notice that we created a random
population with 5 individuals and we had to provide the same number of labels.

Classifier objects are used a callable that receive a single population argument.
It then returns a list of labels corresponding to the assigned classification
of each individual. When we classify the training set, there is a fair chance
of obtaining the original labels:

>>> classifier(pop)                                             # doctest: +SKIP
['A', 'A', 'B', 'B', 'A']

The classifier exposes different classification algorithms that can be accessed
either using the ``pop.classification(labels, <method>)`` method or using the
corresponding attribute ``pop.classification.<method>(labels)``. For instance,
we could try different classifiers

>>> labels = ['A', 'A', 'B', 'B', 'A']
>>> nb = pop.classification.naive_bayes(labels)
>>> svm = pop.classification.svm(labels)

You can check the :cls:`kpop.population.classification.Classification` to see
all available classifiers.


Easy labels
-----------

The default procedure for training a classifier involves passing a list of
labels for the training algorithms. Sometimes those labels can be stored as
meta data in the population object or can be derived from the population
somehow. If the labels argument is a string, kpop will try to obtain the
label list by using the first option valid option:

* Use population.meta[<label>], if it exists.
* If label equals 'ancestry', it creates a list of labels assigning the
id of each sub-population to all its individuals.
* If label is the empty string or None, it looks for a 'labels' column in the
meta information and then returns it.

This interface makes it very convenient to train classifiers to infer population
ancestry. Remember that this is not an admixture analysis since we are assuming
that all individuals belong to a single population.

>>> popA = Population.random(5, 20, id='A')
>>> popB = Population.random(5, 20, id='B')
>>> pop = popA + popB
>>> classifier = pop.classification(labels='ancestry')
>>> classifier(popA)
['A', 'A', 'A', 'A', 'A']
>>> classifier(popB)
['B', 'B', 'B', 'B', 'B']


API docs
--------

.. autoclass:: kpop.population.classification.Classification
    :members:
