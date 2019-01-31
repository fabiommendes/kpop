import re
import sys
from typing import List
from pathlib import Path

STRUCTURE_PARAM_RE = re.compile(r'[#]define +([A-Z0-9]+) +(.+)')


# ----------------------------------------------------------------------------
# Structure data files
#
def read_structure(file, mainparams, extraparams):
    """
    Read data from a Structure project and export to a Kpop population.
    """
    import kpop
    kpop.Population()


def parse_structure(file, mainparams=None, extraparams=None):
    """
    Read data from a Structure and return as dictionary with the following 
    schema::

        { 
            'config': { NAME: VALUE },
            'data': [ {
                'name': NAME,
                'genotype': LOC_DATA, 
            } ]     
        }

    LOC_DATA is a list of lists containing each individual's genotype. 
    """
    file = Path(file)
    mainparams = mainparams or file.parent / 'mainparams'
    extraparams = extraparams or file.parent / 'extraparams'
    config = parse_structure_params(open(mainparams))
    config.update(parse_structure_params(open(extraparams)))

    data = []
    for line in open(file):
        line = line.strip()
        if not line:
            continue

        name, *genotype = line.split()
        genotype = list(map(int, genotype))
        genotype = list(zip(genotype[::2], genotype[1::2]))
        row = {'name': name, 'genotype': genotype}
        data.append(row)

    return {'config': config, 'data': data}


def parse_structure_params(lines: List[str]):
    """
    Read a list of structure configuration lines and return all values as a 
    dictionary.
    """
    result = {}
    for idx, line in enumerate(lines, 1):
        line = line.strip()
        m = STRUCTURE_PARAM_RE.fullmatch(line)
        if m is None and line.startswith('#'):
            msg = 'invalid configuration line: %s) %r' % (idx, line)
            raise ValueError(msg)
        elif m:
            name, value = m.groups()
            result[name] = normalize_structure_param_value(value)
    return result


def normalize_structure_param_value(value):
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


# ----------------------------------------------------------------------------
# PLINK .ped data file
def parse_ped(file):
    """
    Parse a .ped file and return the following data structure::

        [
            {
                'family': FAMILY_ID,
                'id': ID,  # unique within family
                'father_id': ID,
                'mother_id': ID,
                'sex': 'male' | 'female' | None,
                'phenotype': ID,
                'genotype': LOC_DATA,
            }
        ]
    """
    result = []
    for line in open(file):
        family, id_, father, mother, sex, pheno, *data = line.split()
        data = list(map(int, data))
        genotype = list(zip(data[::2], data[1::2]))

        result.append({
            'family': family,
            'id': id_,
            'mother_id': mother,
            'father_id': father,
            'sex': 'male' if sex == '1' else 'female' if sex == '2' else None,
            'phenotype': pheno,
            'genotype': genotype,
        })


# ----------------------------------------------------------------------------
# Converters
def structure_to_ped_data(file, mainparams=None, extraparams=None):
    """
    Convert a Structure project to PLINK's .ped data format.
    """

    result = []
    structure = parse_structure(file, mainparams, extraparams)
    config = structure['config']
    structure_data = structure['data']

    for row in structure_data:
        result.append({
            'family': '0',
            'id': row['name'],
            'mother_id': '0',
            'father_id': '0',
            'sex': None,
            'phenotype': '0',
            'genotype': row['genotype'],
        })
    return result


def structure_to_ped(file, output, mainparams=None, extraparams=None):
    """
    Convert a Structure project to PLINK's .ped + .map files.
    """
    data = structure_to_ped_data(file, mainparams, extraparams)
    save_ped(data, open(output + '.ped', 'w'))
    num_loci = len(data[0]['genotype'])

    with open(output + '.map', 'w') as F:
        for loci in range(num_loci):
            F.write('0 loci%s 0 %s\n' % (loci, loci + 1))


def save_ped(data, file):
    """
    Saves .ped data into file.
    """
    for row in data:
        family = row.get('family', '0')
        id_ = row.get('id', '0')
        father = row.get('father', '0')
        mother = row.get('mother', '0')
        sex = row.get('sex', None)
        sex = '1' if sex == 'male' else '2' if sex == 'female' else '0'
        phenotype = row.get('phenotype', '0')
        genotype = row.get('genotype', [])
        file.write(' '.join([family, id_, father, mother, sex, phenotype]))
        for gen in genotype:
            file.write(' %s %s' % tuple(gen))
        file.write('\n')


if __name__ == '__main__':
    import click
    import json
    from pprint import pprint

    @click.command()
    @click.argument('file')
    @click.option('--from', '-f', help='input file format', prompt=True)
    @click.option('--to', '-t', help='output file format', default='json')
    @click.option('--output', '-o', help='output file path')
    def cli(file, to, output, **kwargs):
        """Perform file conversions"""

        parsers = {'structure': parse_structure, 'ped': parse_ped}
        pairs = {('structure', 'ped'): structure_to_ped}
        from_ = kwargs.get('from', None)

        if to == 'json':
            try:
                parser = parsers[from_]
            except KeyError:
                raise SystemExit('Invalid input file format: %r' % from_)

            data = parser(file)
            json_data = json.dumps(data)
            if not output:
                print(data)
            else:
                with open(output, 'w') as F:
                    F.write(json_data)
        else:
            try:
                converter = pairs[from_, to]
            except KeyError:
                msg = 'invalid combination of files: %r and %r' % (from_, to)
                raise SystemExit(msg)

            result = converter(file, output)
    cli()
