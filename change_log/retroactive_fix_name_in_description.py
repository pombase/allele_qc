import pandas
import glob

manual_changes = pandas.concat([pandas.read_csv(f, sep="\t", na_filter=False) for f in glob.glob("allele_manual*.tsv")])
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
manual_changes.drop(['allele_description'], axis=1, inplace=True)
manual_changes.rename(columns={'change_description_to': 'allele_description'}, inplace=True)
manual_changes[['systematic_id', 'allele_name', 'allele_description', 'change_name_to']].to_csv('retroactive_name_fix.tsv', sep="\t", index=False)