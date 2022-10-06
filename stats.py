# %%
import pandas

data = pandas.read_csv('results/allele_results.tsv', sep='\t', na_filter=False)

aminoacid_data = data[~data['allele_type'].str.contains('nucl')]
fixable_cols = (aminoacid_data['allele_type'] != 'other') & (aminoacid_data['allele_type'] != 'unknown')
fixable_data = aminoacid_data[fixable_cols]
print('fixable_aminoacid_alleles:', len(fixable_data), 1.)
print('needs_fixing:', sum(fixable_data['needs_fixing']), sum(fixable_data['needs_fixing']) / len(fixable_data))

autofix_cols = (fixable_data['rename_to'] != '') | (fixable_data['change_type_to'] != '')
seqerror_cols = fixable_data['sequence_error'] != ''
seqerror_fix_cols = fixable_data['sequence_error'].str.contains(', but')
print('syntax_fixing:', sum(autofix_cols), sum(autofix_cols) / len(fixable_data))
print('sequence_error:', sum(seqerror_cols), sum(seqerror_cols) / len(fixable_data))
print('sequence_error_auto_fix:', sum(seqerror_fix_cols), sum(seqerror_fix_cols) / len(fixable_data))

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
