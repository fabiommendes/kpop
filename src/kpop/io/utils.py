import functools


def file_or_path(func):
    """
    Decorates function so it accepts either file objects or string paths as the
    first argument.
    """

    @functools.wraps(func)
    def decorated(file, *args, **kwargs):
        if isinstance(file, str):
            with open(file) as F:
                result = func(file, *args, **kwargs)
            return result
        else:
            return func(file, *args, **kwargs)

    return decorated
