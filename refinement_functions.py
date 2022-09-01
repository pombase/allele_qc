#%%
import re

def replace_allele_features(regex2modification, pattern_list, matches):
    out_list = list()
    for i in range(len(pattern_list)):
        if type(pattern_list[i]) != str:
            out_list.append(pattern_list[i])
            continue
        if len(matches) == 0:
            for regex_pattern in regex2modification:
                matches += [*re.finditer(regex_pattern, pattern_list[i])]
            matches.sort(key=lambda match : len(match.group()), reverse=True)
        allele_substring = pattern_list[i]
        this_list = [allele_substring]
        for match in matches:
            if match.group() in allele_substring:
                start = allele_substring.find(match.group())
                end = start + len(match.group())
                this_list = [allele_substring[:start], match, allele_substring[end:]]
                # Remove empty strings
                this_list = list(filter(lambda x: x != '', this_list))
                this_list = replace_allele_features(
                    regex2modification, this_list, matches)
                break
        out_list += this_list

    return out_list

def build_regex2modification(modifications):
    out_dict = dict()
    for modification in modifications:
        out_dict[modification['regex']] = modification
    return out_dict

def sort_result(result):
    matches = list()
    unmatched = list()
    for r in result:
        # This is a match
        if type(r) != str:
            matches.append(r)
        # This is a string that has letters or digits in it, so we can't skip it
        elif not re.match('^[^a-zA-Z\d]+$',r):
            unmatched.append(r)
    return matches, unmatched

def allele_is_invalid(allele_description,regex2modification, allele_type, allowed_types):

    result = replace_allele_features(regex2modification, [allele_description], [])
    # Filter out the non-digit non-letter characters
    matches, unmatched = sort_result(result)

    if len(unmatched):
        return 'pattern not matched\t' + ','.join(unmatched)

    expected_list = list()
    invalid_list = list()
    encountered_types = set()
    for match in matches:

        modification = regex2modification[match.re.pattern]
        encountered_types.add(modification['type'])
        invalid_string = modification['check_invalid'](match.groups())
        if invalid_string:
            invalid_list.append(invalid_string)
        else:
            expected_list.append(modification['apply_syntax'](match.groups()))

    if len(invalid_list):
        return 'invalid\t' + ','.join(invalid_list)

    encountered_types = frozenset(encountered_types)
    correct_type = allowed_types[encountered_types]

    corrections = [[],[]]
    if correct_type != allele_type:
        corrections[0].append(f'wrong type')
        corrections[1].append(f'change type to: {correct_type}')

    expected = ','.join(expected_list)
    if expected != allele_description:
        corrections[0].append(f'typo')
        corrections[1].append(f'fix to: {expected}')

    if len(corrections[0]):
        return '\t'.join([','.join(c) for c in corrections])

    return False