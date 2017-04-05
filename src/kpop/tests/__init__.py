import os as _os

# Check if it is running on CI
is_ci = _os.environ.get('CI', 'false') == 'true'