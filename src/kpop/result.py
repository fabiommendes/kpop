import numpy as np
from lazyutils import delegate_ro


#
#  A np.ndarray class that has a __dict__. The constructors are functions that
# create result instances with the proper attributes.
#
class Result(np.ndarray):
    """
    An array sub-type that can hold additional meta data. It is used internally
    by Kpop to store the numerical result of a computation.
    """


def result(data, dtype=None, class_=Result, **kwargs) -> Result:
    """
    Construct a result array from data. All keyword arguments are saved as
    attributes.

    Args:
        data:
            Arbitrary n-dimentional data.

    Returns:
        A result object.
    """

    res = np.asarray(data, dtype=dtype).view(class_)
    for k, v in kwargs.items():
        setattr(res, k, v)
    return res


class Transform(Result):
    """
    Array sub-type that represents the result of an scikit learn transform
    object.
    """

    _data = delegate_ro('transform')

    def transform(self, data):
        """
        Transform the given input data set by the current transformation.
        """

        return self._data.transform(data)


def transform_result(transform: type, data: np.ndarray, **kwargs) -> Transform:
    """
    Apply transform to given data.

    Args:
        transform (type):
            The base class for the transform object (e.g.,
            sklearn.decomposition.PCA)
        data (data):
            Data used to feed the transform.

    Any additional positional and keyword arguments are passed to the transform
    constructor.

    Returns:
        Return a transform array.
    """
    sk_transform = transform(**kwargs)
    res_data = sk_transform.fit_transform(data)
    return result(res_data, class_=Transform, transform=sk_transform)


class Classification(Result):
    pass

def classification_result(classifier, data, **kwargs) -> Classification:
    pass