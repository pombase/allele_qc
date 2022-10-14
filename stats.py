# %%
import pandas

data = pandas.read_csv('results/allele_results.tsv', sep='\t', na_filter=False)

changed_coords = pandas.read_csv('data/gene-coordinate-change-data.tsv', sep='\t', na_filter=False)
changed_coords_ids = set(changed_coords['Systematic ID'])

aminoacid_data = data[~data['allele_type'].str.contains('nucl')]
fixable_cols = (aminoacid_data['allele_type'] != 'other') & (aminoacid_data['allele_type'] != 'unknown')
fixable_data = aminoacid_data[fixable_cols]
print('fixable_aminoacid_alleles:', len(fixable_data), 1.)
print('needs_fixing:', sum(fixable_data['needs_fixing']), sum(fixable_data['needs_fixing']) / len(fixable_data))

autofix_cols = (fixable_data['rename_to'] != '') | (fixable_data['change_type_to'] != '')
seqerror_cols = fixable_data['sequence_error'] != ''
seqerror_fix_cols = fixable_data['sequence_error'].str.contains(', but') & ~fixable_data['systematic_id'].isin(changed_coords_ids)
seqerror_change_coords = seqerror_cols & fixable_data['systematic_id'].isin(changed_coords_ids)
seqerror_left = sum(seqerror_cols) - sum(seqerror_change_coords) - sum(seqerror_fix_cols)

print('syntax_fixing:', sum(autofix_cols), sum(autofix_cols) / len(fixable_data))
print('sequence_error:', sum(seqerror_cols), sum(seqerror_cols) / len(fixable_data))
print('sequence_error_fix_coordinate_change:', sum(seqerror_fix_cols), sum(seqerror_fix_cols) / len(fixable_data))
print('sequence_error_auto_fix:', sum(seqerror_change_coords), sum(seqerror_change_coords) / len(fixable_data))
print('sequence_error_left:', seqerror_left, seqerror_left / len(fixable_data))

print()
nucleotide_data = data[data['allele_type'].str.contains('nucl')]
fixable_cols = (nucleotide_data['allele_type'] != 'other') & (nucleotide_data['allele_type'] != 'unknown')
fixable_data = nucleotide_data[fixable_cols]
print('fixable_nucleotide_alleles:', len(fixable_data), 1.)
print('needs_fixing:', sum(fixable_data['needs_fixing']), sum(fixable_data['needs_fixing']) / len(fixable_data))

autofix_cols = (fixable_data['rename_to'] != '') | (fixable_data['change_type_to'] != '')
seqerror_cols = fixable_data['sequence_error'] != ''
seqerror_fix_cols = fixable_data['sequence_error'].str.contains(', but')
print('syntax_fixing:', sum(autofix_cols), sum(autofix_cols) / len(fixable_data))
print('sequence_error:', sum(seqerror_cols), sum(seqerror_cols) / len(fixable_data))
print('sequence_error_auto_fix:', sum(seqerror_fix_cols), sum(seqerror_fix_cols) / len(fixable_data))


# seqerror_left_log = seqerror_cols & ~(seqerror_fix_cols | seqerror_change_coords)

# for i in fixable_data['allele_description'][seqerror_left_log]:
#     print(i)
