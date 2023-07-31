import pandas
import glob

all_alleles = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
data = pandas.read_csv('../results/allele_cannot_fix_other_errors.tsv', sep='\t', na_filter=False)

already_fixed = pandas.concat([
    pandas.read_csv('manual_cannot_fix_new.tsv', sep='\t', na_filter=False),
    pandas.read_csv('allele_manual_changes_formatted.tsv', sep='\t', na_filter=False),
    pandas.read_csv('cannot_find.tsv', sep='\t', na_filter=False),
    pandas.read_csv('extra_alleles_kim.tsv', sep='\t', na_filter=False),
    ])

# Sometimes the same allele name may have different descriptions
data['unique_id'] = data['allele_name'] + '$' + data['allele_description']
already_fixed['unique_id'] = already_fixed['allele_name'] + '$' + already_fixed['allele_description']

data = data[~data['unique_id'].isin(already_fixed['unique_id'])].copy()

cannot_fix = pandas.read_csv('cannot_find.tsv', sep='\t', na_filter=False)
data = data[~data['allele_name'].isin(cannot_fix['allele_name'])].copy()

# Give a warning for those that are already in the changelog

change_log_files = glob.glob('../change_log/allele*.tsv')
change_log_alleles = pandas.concat([pandas.read_csv(f, sep='\t', na_filter=False) for f in change_log_files])
change_log_alleles['unique_id'] = change_log_alleles['allele_name'] + '$' + change_log_alleles['allele_description']

already_in_changelog = data['unique_id'].isin(change_log_alleles['unique_id'])
warn_already_in_changelog = data[already_in_changelog].copy()
if len(warn_already_in_changelog) > 0:
    print("Warning: the following alleles are already in the changelog:")
    print(warn_already_in_changelog[['allele_name', 'allele_description']])

data = data[~already_in_changelog].copy()

# merge on systematic_id and allele_name
data = pandas.merge(data, all_alleles, on=['systematic_id', 'allele_name', 'allele_description'], how='left')

data.sort_values(['reference', 'systematic_id', 'allele_name', 'allele_description'], inplace=True)

data['change_description_to'] = ''
data['change_name_to'] = ''
data['change_type_to'] = ''
data['comment'] = ''

columns = ['systematic_id', 'allele_name', 'allele_description', 'allele_type', 'reference', 'change_description_to', 'change_name_to', 'change_type_to', 'comment']

data[columns].to_csv('manual_cannot_fix.tsv', sep='\t', index=False)
