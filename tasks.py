import sys
from invoke import run, task
from python_boilerplate.tasks import *


@task
def configure(ctx):
    """
    Instructions for preparing package for development.
    """

    run("{0!s} -m pip install .[dev] -r requirements.txt".format(sys.executable))