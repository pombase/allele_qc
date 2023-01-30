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

# Download data from PomBase and load the genome (last step takes a while)
bash get_data.sh

# Find errors in alleles and propose fixes (see outputs in results/allele_results*.tsv)
python perform_qc.py

## (If correcting coordinates)
# make a list of the alleles that need fixing -> results/alleles_coordinate_change.tsv
python load_coordinate_changes.py
# make a dictionary that contains an alignment of old and new sequences -> results/coordinate_changes_dict.json
python build_alignment_dict.py
# fix the coordinates -> results/allele_results_after_coordinates.tsv
python fix_coordinates.py
## (End of If correcting coordinates)

# Make lists of alleles that can be auto-fixed and those that would require human supervision (results/allele_auto_fix.tsv and results/allele_needs_supervision.tsv)
python get_allele_autofix.py
```

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

To download the data from PomBase, run:

```
bash get_data.sh
```
> **NOTE**: This calls a python script, so remember to also run `poetry shell` before this script. The script uses curl, so if using windows you might have to adapt

This creates the folder `data` and:

* Downloads alleles from PomBase to `data/alleles.tsv` (columns are self-explaining).
* Downloads `.contig` files from latest PomBase release.
* Downloads the `fasta` sequence files from the latest PomBase release.
* Creates a "genome dictionary" (see `load_genome.py`) and stores it in the file `data/genome.pickle`.

## Analysis

### Defining syntax rules

We define "syntax rules" representing the syntax of a type of mutation as dictionaries in a python list that we call a "grammar". The dictionaries are parsed into `SyntaxRule` objects (see [models.py](models.py)).

A full grammar can be found in [grammar.py](grammar.py). Below an example of a rule to represent several single aminoacid mutations, in the form of `VP-120-AA` (Valine and Proline in position 120 and 121 replaced by Alanines).

```python
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

{
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        # This is only valid for cases with two aminoacids or more (not to clash with amino_acid_insertion:usual)
        'regex': f'(?<!\d)({aa}{aa}+)-?(\d+)-?({aa}+)(?!\d)',
        'apply_syntax': lambda g: '-'.join(g).upper(),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else '',
        'check_sequence': lambda g, gg: check_sequence_multiple_pos(g, gg, 'peptide'),
        'coordinate_indexes': (1,)
    },
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

* `coordinate_indexes`: a tuple with the indexes of the regex groups that contain coordinates. This is used to update the allele coordinates if allele sequence changes.

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

> **NOTE**: This requires having successfully ran `get_data.sh`. If no arguments are provided, the analysis scripts read the required files from their default locations, but you can provide different input and output files as arguments (see the docstrings in the scripts).

The main step of the analysis is:

```
python perform_qc.py
```

Takes the allele file as input (by default `data/alleles.tsv`), and generates a new file (by default `results/allele_results.tsv`) that contains the following new columns:

### New columns in allele file

* `allele_parts`: the substrings of the `allele_description` that match regex patterns in the grammar, separated by `|` characters. For example `E325A G338D` would result in `E325A|G338D`.
* `needs_fixing`: `True` or `False` depending on whether the allele needs fixing.
* `change_description_to`: if the correct nomenclature differs from `allele_description`, `change_description_to` contains the right syntax. E.g. for `E325A G338D` contains `E325A,G338D`.
* `rules_applied`: for each of the allele_parts, the syntax rule `type` and `name` as `|`-delimited `type:name`. E.g. for `VP-120-AA,E325A` contains `amino_acid_mutation:multiple_aa|amino_acid_mutation:single_aa`.
* `pattern_error`: contains the parts of the allele that are not picked up by any regular expression in the grammar.
* `invalid_error`: output of the function `check_invalid` of each of the rules applied. E.g. for `KKRKK-71-NEHG` contains `lengths don't match: KKRKK-NEHG`.
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

### Optional - Coordinate changes

Some of the alleles for which `sequence_error`s are found might result from residue coordinates refering to previous gene structures. E.g. if the starting methionine has been changed, all residue coordinates are shifted. To fix this case, we use a genome change log produced with https://github.com/pombase/genome_changelog. Essentially, a tsv file with changes to gene structures. We only do this attempt of fixing for alleles that have mutations in the aminoacid sequence. For DNA sequence, since the probability of getting the right nucleotide by chance is ~25%, we cannot be sure it is safe to switch coordinates even if that gives the right nucleotide.

To make a table of the affected alleles, run:

```
python load_coordinate_changes.py
```

See also the script itself. We consider that alleles of genes for which coordinates were changed can be affected if they are found in publications where sequence errors were found for alleles of the same gene. E.g. for a given gene, a partial deletion allele `30-50` may not give a sequence error by itself. However, if it is in a session with an allele `A58V`, which does give an error, then it will be considered. Then, there are some publications in which only partial deletions exist, in which it may be not possible to tell whether sequence errors exist. For those we cannot know and they are labelled as ambiguous by the script. The script also generates a text file `results/systematic_ids_excluded_coordinate_changes.txt`, for genes of which the sequence coordinates were changed more than once. These will be excluded from the below dictionary.

In order to convert from old coordinates to new coordinates, we build a sequence alignment of both sequences (see docstring). Run:

```
python build_alignment_dict.py
```

Finally, run and see docstring:

```
python fix_coordinates.py
```

### Summarising results

Finally, you can generate two files, one with alleles with errors that can be auto-corrected, and one that may require human input to make the decission.

Run and see docstring:

```
python get_allele_autofix.py
```


## TODO

document mitochondrial pre_fix

## Running the API

```
docker build -t allele_qc_api .
docker run -d --name apicontainer -p 8000:80 allele_qc_api
```