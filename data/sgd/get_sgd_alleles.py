import requests
import sys
import pandas

if __name__ == '__main__':

    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('usage python get_sgd_alleles.py <output_file>')
        sys.exit(0)
    if len(sys.argv) != 2:
        print('Invalid number of arguments')
        sys.exit(1)

    query = """<query model="genomic" view="Gene.secondaryIdentifier Gene.symbol Gene.alleles.name Gene.alleles.aliasName Gene.alleles.description Gene.alleles.alleleClass" sortOrder="Gene.alleles.description ASC" ></query>"""
    # We make the request as json and use pandas because some descriptions contain spaces
    response = requests.get(
        f'https://yeastmine.yeastgenome.org/yeastmine/service/query/results?query={query}&format=json'
    )
    data: pandas.DataFrame = pandas.DataFrame(response.json()['results'])
    data.columns = ['systematic_id', 'primary_name', 'allele_name', 'allele_synonym', 'allele_description', 'allele_type']
    data.allele_description = data.allele_description.str.replace('\s+', ' ', regex=True)

    data.to_csv(sys.argv[1], sep='\t', index=False)
