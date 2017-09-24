import os

from ..libs import np


def csv_lines(data, sep=',', columns=None, align=None,  # noqa: C901
              end=os.linesep, decimal_places=None):
    """
    Renders 2D table ``data`` as CSV. It returns a list of lines that should
    be written in a file object. This function accepts many formatting options
    and a per-column formatting.

    Usage:
        >>> lines = csv_lines([[1, 2], [3, 4]], columns=['A', 'B'])
        >>> print(''.join(lines))
        A,B
        1,2
        3,4

    Args:
        data:
            A uniform sequence of sequences.
        sep:
            The separator between items in the same line.
        columns:
            A list of column names or column specs. If names are not
            given, it will not render the first header column.
        align (str):
            Either None, 'right' or 'left'. Defines the alignment of all
            columns.
        decimal_places (int):
            Number of decimal places used render float data.
        end (str):
            A string appended to the end of each line (can be used to replace
            newlines).

    Returns:
        An iterator over the lines of the generated CSV file.
    """

    # Strategy is to read all columns, format them accordingly and create
    # iterators for each column which are coordinated to create an iterator for
    # each line

    if columns is None:
        columns = [None] * len(data[0])
    cols = list(map(Column, columns))

    # Makes sure that all columns have names if at least one column defines it
    if any(col.name is not None for col in cols):
        for i, col in enumerate(cols, 1):
            col.name = col.name or 'c%s' % i

    # Define alignment
    if align:
        align = 'right' if align is True else align
        for col in cols:
            if col.align is None:
                col.align = align

    # Define the number of decimal places
    if decimal_places:
        for col in cols:
            if col.decimal_places is None:
                col.decimal_places = decimal_places

    # Transpose data
    cols_data = list(np.array(data, dtype=object).T)
    for col, col_data in zip(cols, cols_data):
        col.feed_data(col_data)

    # Create main iterator
    col_iters = list(map(iter, cols))
    while True:
        line = sep.join(map(next, col_iters))
        if not line:
            return
        yield line + end


class Column:
    """
    Represents a column with some formatting options.
    """

    @classmethod
    def from_spec(cls, spec):
        """
        Creates column from spec:
            None -> Column()
            str  -> Column(name=spec)
            dict -> Column(**spec)
        """
        if spec is None:
            return Column()
        elif isinstance(spec, str):
            return Column(name=spec)
        else:
            return Column(**spec)

    def __init__(self, name=None, align=None, decimal_places=None):
        self.name = name
        self.align = align
        self.data = None
        self.decimal_places = decimal_places

    def __iter__(self):
        fmt = self.fmt
        name = None if self.name is None else fmt(self.name)
        data = list(map(fmt, self.data))

        # Reformat data considering align
        if self.align:
            justify = str.ljust if self.align == 'left' else str.rjust
            col_size = max(map(len, data))
            if name is not None:
                col_size = max(col_size, len(name))
                name = justify(name, col_size)
            data = map(lambda x: justify(x, col_size), data)

        if self.name is not None:
            yield name

        for x in data:
            yield x

    def fmt(self, value):
        """
        Format object as an string argument.
        """

        if isinstance(value, float) and self.decimal_places is not None:
            return ('%%.%sf' % self.decimal_places) % value
        return str(value)

    def feed_data(self, data):
        """
        Saves data to be used by the columns iterator.

        Data is only processed when iter(column) is called.
        """

        self.data = data
