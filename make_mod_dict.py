# We could use an ontology library, but we can avoid adding dependencies with a simple script
import re
import json

all_mods_in_pombase = set()

with open('data/pombase-chado.modifications', 'r') as ins:
    for line in ins:
        all_mods_in_pombase.add(line.strip().split('\t')[2])


allowed_mod_dict = dict()

with open('data/mod.obo', 'r') as ins:
    full_file = ins.read()
    regex = re.compile(r'\[Term\](.|\n)+?(?=\n\n)', re.MULTILINE)
    for i in re.finditer(regex, full_file):
        term_dict = dict()
        for element in i.group().split('\n')[1:]:
            key, value = re.match(r'^(\w+): (.+)$', element).groups()
            subvalue_match = re.match(r'^(\w+): "(.+)"$', value)
            if subvalue_match:
                value = subvalue_match.group(2)
                key = subvalue_match.group(1)
            if key not in term_dict:
                term_dict[key] = [value]
            else:
                term_dict[key].append(value)
        # This can be used to skip the mods that are not in annotations already, but we don't do it.
        # if term_dict['id'][0] not in all_mods_in_pombase:
        #     continue

        allowed_mod_dict[term_dict['id'][0]] = ""
        if 'Origin' in term_dict:
            # Keep only aminoacids
            allowed_mod_dict[term_dict['id'][0]] = ''.join(sorted(list(set(x for x in term_dict['Origin'][0].split(', ') if x in 'GPAVLIMCFYWHKRQNEDST'))))

with open('data/allowed_mod_dict.json', 'w') as outfile:
    json.dump(allowed_mod_dict, outfile, indent=4, sort_keys=True)
