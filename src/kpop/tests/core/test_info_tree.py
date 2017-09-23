import pytest

from kpop.info_tree import InfoTree


class TestInfoTree:
    """
    Test the info tree data structure.
    """

    @pytest.fixture
    def tree(self):
        return InfoTree([('foo', 42), ('bar', 0)])

    @pytest.fixture
    def nested_tree(self, tree):
        tree.bar = InfoTree([('bar', 1)])
        return tree

    def test_info_tree_has_basic_attr_interface(self, tree):
        assert tree.foo == 42
        assert tree.bar == 0
        tree.baz = 1
        assert tree.baz == 1

    def test_info_tree_repr(self, tree):
        assert repr(tree) == (
            'foo: 42\n'
            'bar: 0'
        )

    def test_nested_tree(self, nested_tree):
        assert repr(nested_tree) == (
            'foo: 42\n'
            'bar:\n'
            '  bar: 1'
        )

    def test_to_dict(self, nested_tree):
        assert nested_tree.to_dict() == {'foo': 42, 'bar': {'bar': 1}}

    def test_to_json(self, nested_tree):
        print(nested_tree.to_json())
        assert nested_tree.to_json() == (
            '{\n'
            '  "foo": 42,\n'
            '  "bar": {\n'
            '    "bar": 1\n'
            '  }\n'
            '}'
        )

    def test_info_tree_errors(self, tree):
        with pytest.raises(AttributeError):
            tree.baz

        with pytest.raises(AttributeError):
            del tree.baz
