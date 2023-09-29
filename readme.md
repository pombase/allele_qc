# Allele QC for PomBase

A series of scripts for quality control of alleles. Currently used for PomBase, but could be adapted.

It checks that:

* Allele descriptions adhere to the allele nomenclature.
* Coordinates used in the allele descriptions are within the sequence length (E.g. we can't say a protein is missing aminoacids 120-150 if the protein length is 140).
* Aminoacids or nucleotides mentioned in allele descriptions occur at the indicated positions (E.g. if we say alanin 120 is mutated to glycine, check that the aminoacid at position 120 is indeed alanine).

## TL;DR;

```bash
# Install dependencies
poetry install

# Activate python environment
poetry shell

# Set up the necessary transvar variables (you must have installed transvar, see next section)
. transvar_env_vars.sh
bash set_up_transvar.sh


# Run this script (See the comments in the subscripts)
bash run_analysis.sh
```

## Installing

### Python dependencies

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

### Other dependencies and setting up transvar

This project uses [transvar](https://github.com/zwdzwd/transvar). This requires to install some binaries.

```bash

# If you have linux and you want to install them globally
sudo apt install -y samtools tabix

# If you want to install them locally (see the content of the script)
# > basically downloads the libs and uses make to build the necessary bin files, then deletes all unnecesary source code
bash install_transvar_dependencies_locally.sh

```

Then, regardless of whether you are using local or global installation of `samtools` and `tabix`:

```bash
# Env vars (see script)
. transvar_env_vars.sh

# Build the transvar database, and test that it works
bash set_up_transvar.sh
```

## What the pipeline does

The best thing is to look at the script `run_analysis.sh`, the subscripts are well documented.

### Defining syntax rules in a grammar

These are used to interpret the allele descriptions, check that the sequence residues they refer to are correct, and to format the description correctly/

We define "syntax rules" representing the syntax of a type of mutation as dictionaries in a python list that we call a "grammar". The dictionaries are parsed into `SyntaxRule` objects (see [models.py](models.py)).

A full grammar can be found in [grammar.py](grammar.py), and the best is to go through that example and the tests to understand how it works. Below an example of a rule to represent several single aminoacid mutations, in the form of `VP120AA` (Valine and Proline in position 120 and 121 replaced by Alanines).

```python
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

{
        'type': 'amino_acid_mutation',
        'rule_name': 'single_aa',
        'regex': f'(?<=\\b)({aa})(\d+)({aa})(?=\\b)',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'further_check': lambda g, gg: g[0] != g[2],
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}{g[1]}{g[2]}']
    },
```

* `type`: the type of mutation (see below)
* `rule_name`: `type` + `rule_name` should be unique, it identifies the rule matching a particular mutation.
* `regex`: a regular expression with a pattern that represents a permisive syntax pattern for this type of mutation.
  * Note the use of `{aa}` inside the python f string to represent any aminoacid.
  * The pattern is a regex, so it can be permissive
  * The capture groups of the regex are passed to the functions below as a first argument.
* `apply_syntax`: a function that takes a tuple of `re.Match[str]`, representing a match in a string to the pattern described in `regex`, and returns the correctly-formatted mutation.
  * In the example `VP120AA`, the groups are `('VP', '120', 'AA')`, and the function returns `VP120AA`.
  * The function can be defined inline using `lambda`, or outside of the dictionary.
* `check_sequence`: a function that takes the two arguments below and checks whether the sequence position indicates exist (the index is within the boundaries of the sequence), and contain the indicated residue. If sequence is correct, should return empty string, otherwise, the residue or residues that are incorrect, see the examples in `grammar.py` and the tests.
  * A tuple of `re.Match[str]`, representing a match in a string to the pattern described in `regex`
  * A "gene dictionary", see `load_genome.py`
* `further_check`: takes the same argument as above, and is used as an extra check on top of matching the `regex` pattern. In this example, it checks that the parts before and after the number are different (e.g. `VP120VP` will not be matched by this grammar rule, even if it matches the pattern).
* `format_for_transvar`: takes the same arguments as above, and returns the transvar formatted variant.


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

### Config file

See `config.json`, the variables included there are self-explaining. See also `data/sgd/config.sgd.json` for SGD.

### New columns in allele file after analysis

* `allele_parts`: the substrings of the `allele_description` that match regex patterns in the grammar, separated by `|` characters. For example `E325A G338D` would result in `E325A|G338D`.
* `needs_fixing`: `True` or `False` depending on whether the allele needs fixing.
* `change_description_to`: if the correct nomenclature differs from `allele_description`, `change_description_to` contains the right syntax. E.g. for `E325A G338D` contains `E325A,G338D`.
* `rules_applied`: for each of the allele_parts, the syntax rule `type` and `name` as `|`-delimited `type:name`. E.g. for `VP120AA,E325A` contains `amino_acid_mutation:multiple_aa|amino_acid_mutation:single_aa`.
* `pattern_error`: contains the parts of the allele that are not picked up by any regular expression in the grammar.
* `invalid_error`: if the systematic_id or sequence of a gene product is missing.
* `sequence_error`: contains an error if the indicated position does not exist or if the aminoacid indicated at a certain position is not correct. Values to be homogenised in the future.
* `change_type_to`: if `allele_type` is not right, contains the right type. This comes from using the `frozenset` in `grammar.py`.
* `auto_fix_comment`: for alleles that can be auto-fixed, contains some info (e.g.):
  * `syntax_error`
  * `type_error`
  * `syntax_and_type_error`
  * `multi_shift_fix`: For a given publication with 4 or more alleles with errors, if shifting all coordinates by the same amount fixes the error, the fix is accepted. This fix has the lowest priority.
  * `histone_fix`: histone indexes have been often counted ignoring the first methionine, so if a shift by +1 in index fixes the sequence error in a histone, it is accepted. This type of fix takes second priority.
  * `old_coords_fix, revision xxx: gene_coordinates`: if using old gene coordinates the error is fixed, the fix is accepted. This type of fix takes max priority, as it is in principle more reliable.
* `solution_index`: Normally empty, but if more than one solution has been found to a sequence error, they have the index of the solution.

### Optional - Using old coordinta changes for fixes

Some of the alleles for which `sequence_error`s are found might result from residue coordinates refering to previous gene structures. E.g. if the starting methionine has been changed, all residue coordinates are shifted. To fix this case, we use a genome change log produced with https://github.com/pombase/genome_changelog (for PomBase, see  `build_alignment_dict_from_genome.py`). For DNA sequence, since the probability of getting the right nucleotide by chance is ~25%, we cannot be sure it is safe to switch coordinates even if that gives the right nucleotide.

For sgd data, we use a method that uses only the protein sequences: https://github.com/pombase/all_previous_sgd_peptide_sequences, see `build_alignment_dict_from_peptides.py`.

## Running the API in Docker

```
docker build -t allele_qc_api .
docker run -d --name apicontainer -p 8000:80 allele_qc_api
```

Then if you go to http://localhost:8000/ you should be redirected to the API documentation and
you can run a test request directly there.
