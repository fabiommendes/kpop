class Attr:
    """
    Implements special attributes.
    """

    def __init__(self, pop: 'kpop.Population' = None):
        if isinstance(pop, type):
            self._population = None
            self._cls = pop
        else:
            self._population = pop
            self._cls = type(pop)

    def __get__(self, instance, cls=None):
        if instance is None:
            return type(self)(cls)
        else:
            attr = type(self)(instance)
            setattr(instance, self.__class__.__name__.lower(), attr)
            return attr

    def _new_population(self, *args, **kwargs):
        from kpop import Population
        return Population(*args, **kwargs)