import os

from jupyter_client import KernelManager
from notebook.notebookapp import NotebookApp

DIRNAME = os.path.dirname(__file__)


class KpopKernelManager(KernelManager):
    kernel_name = 'python'
    kernel_script_path = os.path.join(DIRNAME, 'jupyter_kernel.py')


class ZMQTerminalKpopApp(NotebookApp):
    name = 'kpop'
