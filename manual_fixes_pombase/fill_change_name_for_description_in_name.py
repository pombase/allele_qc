import pandas
import sys

manual_changes = pandas.read_csv(sys.argv[1], sep="\t", na_filter=False)
original_dataset = manual_changes.copy()
column_order = list(original_dataset.columns)
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
print(original_dataset.shape)
final_dataset = original_dataset.merge(manual_changes[['allele_name', 'change_name_to']], on='allele_name', how='left')
final_dataset.fillna('', inplace=True)
final_dataset['change_name_to'] = final_dataset.apply(lambda row: row['change_name_to_y'] if row['change_name_to_y'] != '' else row['change_name_to_x'], axis=1)
final_dataset.drop(['change_name_to_x', 'change_name_to_y'], axis=1, inplace=True)
final_dataset[column_order].to_csv('names_filled.tsv', sep="\t", index=False)