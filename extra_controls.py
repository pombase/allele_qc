import pandas
import re
data = pandas.read_csv('results/allele_auto_fix.tsv', delimiter='\t', na_filter=False)

print('\033[0;32mmixed case\033[0m')
print()
for i, row in data.iterrows():
    desc = row['allele_description']
    # Combinations of upper and lower case might have been weird patterns
    if re.search('[a-z]', desc) and re.search('[A-Z]', desc):
        print(desc)

print()
print('\033[0;32mnt or aa\033[0m')
print()
for i, row in data.iterrows():
    desc = row['allele_description']

    # Presence of aa or nt might mean nucleotides or aminoacids
    if 'aa' in desc.lower() or 'nt' in desc.lower():
        print(desc)
