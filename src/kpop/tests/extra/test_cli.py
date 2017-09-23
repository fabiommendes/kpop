import io
import sys

import pytest

from kpop.cli import main

pytestmark = [pytest.mark.slow]


def exec_main(argv):
    try:
        old_out, sys.stdout = sys.stdout, io.StringIO()
        old_argv = sys.argv
        sys.argv = argv
        main()
    except SystemExit:
        pass
    finally:
        data = sys.stdout.getvalue()
        sys.argv = old_argv
        sys.stdout = old_out
    return data


class TestImport:
    def test_import_command(self, path):
        for file in ['popAB.csv', 'popAB.pickle']:
            path_ = path(file)
            out = exec_main(['kpop', 'import', path_]).strip()
            print(out)
            assert out.startswith('Population data')
            assert 'A1: 11 22 12 12 12' in out
            assert 'B1: 22 22 22 12 22' in out
