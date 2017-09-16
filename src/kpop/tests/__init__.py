import os as _os

# Check if it is running on CI
is_travis = _os.environ.get('TRAVIS', 'false') != 'false'
is_appveyor = _os.environ.get('APPVEYOR', 'false') != 'false'
is_ci = is_travis or is_appveyor or _os.environ.get('CI', 'false') != 'false'

def load_data(name):
    return open(_os.path.join(_os.path.dirname(__file__), 'data', name))