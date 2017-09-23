import os as _os
from contextlib import contextmanager as _ctxmanager

# Check if it is running on CI
is_travis = _os.environ.get('TRAVIS', 'false') != 'false'
is_appveyor = _os.environ.get('APPVEYOR', 'false') != 'false'
is_ci = is_travis or is_appveyor or _os.environ.get('CI', 'false') != 'false'


def load_data(name):
    return open(data_path(name))


def data_path(name):
    return _os.path.join(_os.path.dirname(__file__), 'data', name)


@_ctxmanager
def temporary_location():
    curr_path = _os.getcwd()
    tmp_path = _os.path.join(_os.path.dirname(__file__), 'tmp')
    _os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        _os.chdir(curr_path)
