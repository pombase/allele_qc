#%%
import re

def replace_allele_features(syntax_rules, pattern_list, matches):
    out_list = list()
    for i in range(len(pattern_list)):
        if type(pattern_list[i]) != str:
            out_list.append(pattern_list[i])
            continue
        if len(matches) == 0:
            for syntax_rule in syntax_rules:
                matches += [*re.finditer(syntax_rule['regex'], pattern_list[i])]
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
                    syntax_rules, this_list, matches)
                break
        out_list += this_list

    return out_list

def build_regex2syntax_rule(syntax_rules):
    """
    A dictionary in which the
    """
    out_dict = dict()
    for syntax_rule in syntax_rules:
        out_dict[syntax_rule['regex']] = syntax_rule
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

def allele_is_invalid(allele_description,syntax_rules, allele_type, allowed_types, gene):

    result = replace_allele_features(syntax_rules, [allele_description], [])
    # Filter out the non-digit non-letter characters
    matches, unmatched = sort_result(result)

    if len(unmatched):
        return 'pattern not matched\t' + ','.join(unmatched) + '\t'

    regex2syntax_rule = build_regex2syntax_rule(syntax_rules)
    expected_list = list()
    invalid_list = list()
    sequence_error_list = list()
    encountered_types = set()
    applied_rules = list()
    for match in matches:

        syntax_rule = regex2syntax_rule[match.re.pattern]
        encountered_types.add(syntax_rule['type'])
        invalid_string = syntax_rule['check_invalid'](match.groups())
        if invalid_string:
            invalid_list.append(invalid_string)
            continue

        sequence_error = syntax_rule['check_sequence'](match.groups(), gene)
        if sequence_error:
            sequence_error_list.append(sequence_error)

        expected_list.append(syntax_rule['apply_syntax'](match.groups()))
        applied_rules.append(f'{syntax_rule["type"]}:{syntax_rule["rule_name"]}')

    error_list = list()
    if len(invalid_list):
        error_list.append(','.join(invalid_list))
    if len(sequence_error_list):
        error_list.append(','.join(sequence_error_list))
    if len(error_list):
        return 'invalid\t' + '|'.join(error_list)

    encountered_types = frozenset(encountered_types)
    correct_type = allowed_types[encountered_types]

    corrections = [[],[]]
    if correct_type != allele_type:
        corrections[0].append(f'wrong type')
        corrections[1].append(f'change type to: {correct_type}')

    expected = ','.join(expected_list)
    if expected != allele_description:
        corrections[0].append(f'typo:')
        corrections[1].append(f'fix to: {expected}')

    corrections.append(applied_rules)
    if len(corrections[0]):
        return '\t'.join([','.join(c) for c in corrections])

    return False