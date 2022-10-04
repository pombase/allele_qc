# Allele QC for PomBase

A series of scripts for quality control of alleles in PomBase. It checks that:

* Allele descriptions adhere to the allele nomenclature.
* Coordinates used in the allele descriptions are within the sequence length (E.g. we can't say a protein is missing aminoacids 120-150 if the protein length is 140).
* Aminoacids or nucleotides mentioned in allele descriptions occur at the indicated positions (E.g. if we say alanin 120 is mutated to glycine, check that the aminoacid at position 120 is indeed alanine).

## Installing

To install the dependencies, we used poetry (see [poetry installation instructions](https://python-poetry.org/docs/)).

In the source directory run:

```
poetry install
```

This should create a folder `.venv` with the python virtual environment. To activate the virtual environment, then run:

```
poetry shell
```

Now when you call `python`, it will be the one from the `.venv`.

## Getting the data

Next, we need to download the data from PomBase, for that run (this script uses curl, so if using windows you might have to adapt):

```
bash get_data.sh
```
> **NOTE**: This calls a python script, so remember to also run `poetry shell` before this script.

This creates the folder `data` and:

* Downloads alleles from PomBase to `data/alleles.tsv` (columns are self-explaining).
* Downloads `.contig` files from latest PomBase release.
* Creates a "genome dictionary" (see `load_genome.py`) and stores it in the file `data/genome.pickle`.

## Allele syntax

TODO, mention also chained mutations, comma-separated.

## Analysis

### Defining syntax rules

We define "syntax rules" representing the syntax of a type of mutation as dictionaries in a python list that we call a "grammar" (see `grammar.py`). Below an example of a rule to represent single aminoacid mutations, in the form of `VP-120-AA` (Valine and Proline in position 120 and 121 replaced by Alanines).

```python
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

{
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        # This is only valid for cases with two aminoacids or more (not to clash with amino_acid_insertion:usual)
        'regex': f'({aa}{aa}+)-?(\d+)-?({aa}+)(?!\d)',
        # We fix the case in which dashes are used for a single aa substitution: K-90-R
        'apply_syntax': lambda g: '-'.join(g).upper() if len(g[0])!=1 else ''.join(g).upper(),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else False,
        'check_sequence': check_sequence_multiple_aa
}
```

* `type`: the type of mutation (see below)
* `rule_name`: `type` + `rule_name` should be unique, it identifies the rule matching a particular mutation.
* `regex`: a regular expression with a pattern that represents a permisive syntax pattern for this type of mutation.
  * Note the use of `{aa}` inside the python f string to represent any aminoacid.
  * The pattern is permissive, because it will match the right syntax `VP-120-AA`, but also `VP120AA`)
* `apply_syntax`: a function that takes a tuple of `re.Match[str]`, representing a match in a string to the pattern described in `regex`, and returns the correctly-formatted mutation.
  * In the example `VP120AA`, the groups are `('VP', '120', 'AA')`, and the function returns `VP-120-AA`.
  * The function can be defined inline using `lambda`, or outside of the dictionary.
* `check_invalid`: a function that takes a tuple of `re.Match[str]`, representing a match in a string to the pattern described in `regex`, and returns a string describing an error, if a formatting error exists.
  * In the example of `amino_acid_mutation:multiple_aa` an error is returned if the number of aminoacids before and after the number does not match, e.g. `CP120AAA`.
* `check_sequence`: a function that takes two arguments:
  * A tuple of `re.Match[str]`, representing a match in a string to the pattern described in `regex`
  * A "gene dictionary", see `load_genome.py`

  The function verifies that the proposed mutation is compatible with gene DNA or peptide sequence, and returns an error string otherwise (see examples in `grammar.py`).

### Defining allele categories

In PomBase we use categories for allele types, depending on the types of mutations they contain. These are described in a dictionary that uses `frozenset` objects as keys in `grammar.py`. For example:

```python
allowed_types = {
    frozenset({'amino_acid_mutation'}): 'amino_acid_mutation',
    frozenset({'partial_amino_acid_deletion'}): 'partial_amino_acid_deletion',
    frozenset({'amino_acid_mutation','partial_amino_acid_deletion'}): 'partial_amino_acid_deletion'
    }
```

This is convenient because it allows to represent the fact that an allele that contains only the types of mutations indicated by a frozenset in the dictionary is of that type, e.g., if it has only `amino_acid_mutation`, it is of type `amino_acid_mutation`, if it contains `amino_acid_mutation` and `partial_amino_acid_deletion`, it is of type `amino_acid_deletion_and_mutation`.

### Running the analysis

You can run the analysis by:

> **NOTE**: This requires having successfully ran `get_data.sh`

```
python perform_qc.py
```

See `perform_qc.py` 


