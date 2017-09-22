import subprocess


def whereis(bin):
    """
    Returns the path to the executable or None if executable is not found.
    """

    return None


def run_or_quit(cmd, *args, message=None, job_dir=None):
    """
    Run a command in shell with the given args.

    If the command does not exists, raises a RuntimeError that displays the
    given message.

    Args:
        message (str):
            An error message displayed when command is not found.

    Usage:
        >>> run_or_quit('plink', '--noweb', '--file', 'test')
    """

    bin_path = whereis(cmd)

    if bin_path is None:
        if not message:
            message = (
                '"%s" was not found in your system or it does not seem to be '
                'present in the $PATH. Please install it and make it available '
                'in the system\' path.' % cmd
            )
        raise RuntimeError(message)

    error = subprocess.check_call(cmd, shell=True, cwd=job_dir)
    if error != 0:
        raise RuntimeError('structure returned with error code: %s' % error)
