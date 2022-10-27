import pandas

data = pandas.read_csv('results/allele_auto_fix.tsv', sep='\t', na_filter=False)

print()
print('\033[0;32m', '> auto-fix alleles: ', '\033[0m', len(data), sep='')

for f in pandas.unique(data.fix_type):
    print('  ', f, sum(data.fix_type == f))

data = pandas.read_csv('results/allele_needs_supervision.tsv', sep='\t', na_filter=False)

print()
print('\033[0;32m', '> need_supervision: ', '\033[0m', len(data), sep='')

for f in pandas.unique(data.fix_type):
    print('  ', f, sum(data.fix_type == f))

data = pandas.read_csv('results/allele_cannot_fix.tsv', sep='\t', na_filter=False)

print()

nt_number = sum(data['allele_type'].str.contains('nucl'))
aa_number = sum(data['allele_type'].str.contains('amino'))
other_number = sum(~data['allele_type'].str.contains('amino') & ~data['allele_type'].str.contains('nucl'))

invalid_error = sum(data['invalid_error'] != '')
sequence_error = sum(data['sequence_error'] != '')
pattern_error = sum(data['pattern_error'] != '')

print('\033[0;32m', '> cannot fix alleles: ', '\033[0m', len(data), sep='')
print('   nucleotide alleles:', nt_number)
print('   aminoacid alleles:', aa_number)
print('   disruption alleles:', other_number)
print()
print('   sequence_error:', sequence_error, 'of which', sum((data['sequence_error'] != '') & data['allele_type'].str.contains('nucl')), 'nucleotide')
print('   pattern_error:', pattern_error)
print('   invalid_error:', invalid_error)

print()
