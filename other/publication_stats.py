import pandas
from matplotlib import pyplot as plt
import numpy as np

data = pandas.read_csv('../data/phenotype_annotations.phaf', sep='\t', na_filter=False).loc[:, ['Gene systematic ID', 'Allele name', 'Reference']]
# rename columns
data.columns = ['systematic_id', 'allele_name', 'reference']

# data.loc[:, 'reference'] = data['reference'].apply(str.split, args=[','])
# data = data.explode('reference')

# Get all unique pairs of reference and systematic_id
genes = data[['reference', 'systematic_id']].drop_duplicates()
alleles = data[['reference', 'allele_name']].drop_duplicates()

# Count the number of occurrences in each reference
gene_counts = genes.groupby('reference').count()
allele_counts = alleles.groupby('reference').count()

# Exclude genome-wide studies
gene_counts = gene_counts[gene_counts['systematic_id'] < 100]
gene_counts['reference'] = gene_counts.index
gene_counts['pmid'] = gene_counts['reference'].apply(lambda x: int(x.split(':')[1]))

allele_counts = allele_counts[allele_counts['allele_name'] < 100]
allele_counts['reference'] = allele_counts.index
allele_counts['pmid'] = allele_counts['reference'].apply(lambda x: int(x.split(':')[1]))

print(gene_counts)

plt.figure()
plt.hist(gene_counts['systematic_id'], bins=range(0, 30))
plt.title(f'genes per publication, mean {np.mean(gene_counts["systematic_id"])}, median {np.median(gene_counts["systematic_id"])}')
plt.ylabel('Number of publications')
plt.xlabel('Alleles per publication')


plt.figure()
plt.hist(allele_counts['allele_name'], bins=range(0, 30))
plt.title(f'alleles per publication, mean {np.mean(allele_counts["allele_name"])}, median {np.median(allele_counts["allele_name"])}')
plt.ylabel('Number of publications')
plt.xlabel('Alleles per publication')

plt.figure()
plt.scatter(gene_counts['pmid'], gene_counts['systematic_id'])
# plot the regression line
z = np.polyfit(gene_counts['pmid'], gene_counts['systematic_id'], 1)
p = np.poly1d(z)
plt.plot(gene_counts['pmid'], p(gene_counts['pmid']), "r--")


plt.figure()
plt.scatter(allele_counts['pmid'], allele_counts['allele_name'])
# plot the regression line
z = np.polyfit(allele_counts['pmid'], allele_counts['allele_name'], 1)
p = np.poly1d(z)
plt.plot(allele_counts['pmid'],p(allele_counts['pmid']),"r--")

plt.show()
