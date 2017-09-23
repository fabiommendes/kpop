from ..libs import sk_base


def is_sklearn_transformer(obj):
    """
    Return True if obj has a basic Scikit Learn transform API.
    """
    return isinstance(obj, type) and hasattr(obj, 'fit_transform')


def is_sklearn_classifier(obj):
    """
    Return True if obj has a basic Scikit Learn classifier API.
    """
    return sk_base.is_classifier(obj)
