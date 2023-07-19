import pandas
import sys

manual_changes = pandas.read_csv(sys.argv[1], sep="\t", na_filter=False)
manual_changes.fillna('', inplace=True)

# only the ones that did not have a value in change_name_to
manual_changes = manual_changes[(manual_changes['change_name_to'] == '') & (manual_changes['change_description_to'] != '') & (manual_changes['allele_description'] != '')].copy()


def get_change_name_to(row):
    new_description = row['change_description_to']
    old_description = row['allele_description']
    name = row['allele_name']
    if old_description in name:
        return name.replace(old_description, new_description)
    return ''


manual_changes['change_name_to'] = manual_changes.apply(get_change_name_to, axis=1)
manual_changes = manual_changes[manual_changes['change_name_to'] != ''].copy()
manual_changes.to_csv('names_that_should_be_fixed.tsv', sep="\t", index=False)