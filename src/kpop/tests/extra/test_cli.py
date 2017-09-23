import io
import os
import sys

from mock import patch

from kpop import Population
from kpop.cli import main
from kpop.population.io import Io


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


def save_to_list(lst):
    def save(obj, path, format='auto', **kwargs):
        kwargs['format'] = format
        lst.append((obj._population, path, kwargs))

    return patch.object(Io, 'save', save)


class TestImport:
    def test_create_command(self, path):
        lst = []
        with save_to_list(lst):
            exec_main(
                ['kpop', 'create', '-s', '5', '-l', '8', '-i', 'A', 'file.csv'])
        pop, path, kwargs = lst[0]
        assert pop.shape == (5, 8, 2)
        assert pop.id == 'A'
        assert path == 'file.csv'
        assert isinstance(pop, Population)

    def test_export_command(self, path, temp_path, popA):
        with temp_path() as tmp:
            os.unlink(os.path.join(tmp, 'popAB', 'data.txt'))
            os.unlink(os.path.join(tmp, 'popAB', 'mainparams'))
            os.unlink(os.path.join(tmp, 'popAB', 'extraparams'))

            file = path('popAB.pickle')
            cmd = ['kpop', 'export', '-f', 'structure', '-o', 'popAB', file]
            exec_main(cmd)

            # TODO: fix the label
            with open(os.path.join(tmp, 'popAB', 'data.txt')) as F:
                assert F.readline() == 'A1  1 1 2 2 1 2 1 2 1 2\n'

    def test_import_command(self, path):
        for file in ['popAB.csv', 'popAB.pickle']:
            path_ = path(file)
            out = exec_main(['kpop', 'import', path_]).strip()
            print(out)
            assert out.startswith('Population data')
            assert 'A1: 11 22 12 12 12' in out
            assert 'B1: 22 22 22 12 22' in out

    def test_import_and_save_command(self, path, temp_path):
        with temp_path() as tmp:
            path_ = path('popAB.pickle')
            out = exec_main(['kpop', 'import', path_, '-o', 'test.csv']).strip()
            with open('test.csv') as dest, open(path('popAB.csv')) as src:
                for line_dest, line_src in zip(dest, src):
                    assert line_dest == line_src
