import sys
from ipykernel.ipkernel import IPythonKernel
from ipykernel.kernelapp import IPKernelApp
from ipykernel.zmqshell import ZMQInteractiveShell
from lazyutils import lazy
from traitlets import Type

import kpop


class KpopShell(ZMQInteractiveShell):
    def init_user_ns(self):
        """
        Additional symbols for the shell environment.
        """

        super().init_user_ns()
        ns = self.user_ns


class KpopKernel(IPythonKernel):
    implementation = 'kpop'
    implementation_version = kpop.__version__
    language = 'Python 3.x (Kpop)'
    language_version = sys.version

    @lazy
    def banner(self):
        return self.transpyler.get_console_banner()

    @lazy
    def language_info(self):
        transpyler = self.transpyler
        return {
            'mimetype': 'application/x-python',
            'file_extension': '.py',
            'codemirror_mode': {
                "version": 3,
                "name": "kpop"
            },
            'pygments_lexer': 'python',
        }

    shell_class = Type(KpopShell)


def start_kernel():
    """
    Start Pytuga Jupyter kernel.
    """

    IPKernelApp.launch_instance(kernel_class=KpopKernel)


if __name__ == '__main__':
    start_kernel()
