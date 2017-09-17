def is_transformer(obj):
    """
    Return True if obj respects a basic Scikit Learn transform API.
    """
    return isinstance(obj, type) and hasattr(obj, 'fit_transform')