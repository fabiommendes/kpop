class DescriptorMeta(type):
    """
    Metaclass for attribute descriptors.
    """
    def __get__(self, instance, cls=None):
        if instance is None:
            return self
        else:
            return self(instance)