# Then you can do cat mods_with_names.tsv| cut -f3,10|grep phosph|sort|uniq > phospho.txt

import re
import pandas
all_mods_in_pombase = set()

data = pandas.read_csv('../data/pombase-chado.modifications', sep='\t', na_filter=False)
data.columns = ['systematic_id', 'primary_name', 'modification', 'evidence', 'sequence_position', 'annotation_extension', 'reference', 'taxon', 'date']

term_dict = dict()
with open('../data/mod.obo', 'r') as ins:
    full_file = ins.read()
    regex = re.compile(r'id: (.+?)\nname: (.+?)\n', re.MULTILINE)
    for i in re.finditer(regex, full_file):
        id, name = i.groups()
        term_dict[id.strip()] = name.strip()

data['mod_name'] = data['modification'].apply(lambda x: term_dict[x])

data.to_csv('../mods_with_names.tsv', sep='\t', index=False)