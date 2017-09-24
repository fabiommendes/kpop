import json
from collections import OrderedDict
from os import linesep


class InfoTree:
    """
    InfoTree objects are a convenient way to organize information easily
    reportable to the user.
    """
    __slots__ = '_data'

    def __init__(self, data=()):
        super().__setattr__('_data', OrderedDict(data))

    def __repr__(self):
        lines = []
        for k, v in self._data.items():
            if isinstance(v, InfoTree):
                lines.append('%s:' % k)
                lines.append(indent(repr(v), '  '))
            else:
                lines.append('%s: %s' % (k, v))
        return linesep.join(lines)

    def __getattr__(self, item):
        try:
            return self._data[item]
        except KeyError:
            raise AttributeError(item)

    def __setattr__(self, item, value):
        self._data[item] = value

    def __delattr__(self, item):
        try:
            del self._data[item]
        except KeyError:
            raise AttributeError(item)

    def __eq__(self, other):
        if isinstance(other, InfoTree):
            return self._data == other._data
        return NotImplemented

    def to_dict(self, dict=dict):
        """
        Convert tree to a nested dictionary structure.

        You can pass a different dictionary-like class to be used instead of
        dict.
        """
        out = dict(self._data)
        for k, v in out.items():
            if isinstance(v, InfoTree):
                out[k] = v.to_dict()
        return out

    def to_json(self):
        """
        Renders as a JSON string.
        """
        data = self.to_dict(OrderedDict)
        return json.dumps(data, indent=2)

    def to_yaml(self):
        """
        Renders as an YAML string.
        """
        raise NotImplementedError

    def pprint(self):
        """
        Prints to the terminal using pretty print hints such as bold faces and
        colors.
        """
        raise NotImplementedError


def indent(st, indent, sep=linesep):
    """
    Prepend the indent to all lines in the given string.
    """
    return sep.join(indent + line for line in st.split(sep))
