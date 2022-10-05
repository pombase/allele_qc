# %%
import re


def replace_substring_by_match(input_str: str, match: re.Match) -> list[str, re.Match]:
    """
    Convert string into splitted list using match. E.g.:

    input_str = 'A321B'
    match = <re.Match for '321'>
    returns: ['A', <re.Match for '321'>, 'B']
    """
    start = input_str.find(match.group())
    end = start + len(match.group())
    this_list = [input_str[:start], match, input_str[end:]]
    # Remove empty strings
    return list(filter(lambda x: x != '', this_list))


def replace_allele_features(regex_patterns: list[str], input_list: list[str, re.Match], matches: list[re.Match]) -> list[str, re.Match]:
    """
    Looks for matches to the regex patterns in `regex_patterns` in the strings in `input_list`,
    if `matches` is an empty list. If `matches` is not empty, it uses those matches.

    Then, for each match, starting from the longest one, it splits the strings of `input_list` into substrings
    and a match object. For example, for regex: \d+ applied to `input_list` ['V320A'], it would return ['V', Match Object matching 320, 'A'].

    The function is recursive, since `input_list` changes every time that a match is substituted.

    Example input:

    regex_patterns = ['\d+', '[a-zA-Z]']
    input_list = ['A321B**']
    matches = []

    returns: [<re.Match for 'A'>, <re.Match for '321'>, <re.Match for 'B'>, '**']
    """
    # The output, that will be identical to input_list if no pattern is found.
    out_list = list()

    for allele_substring in input_list:

        # If the element is a re.Match, we include it as is.
        if allele_substring != str:
            out_list.append(allele_substring)
            continue

        # If matches are not provided, we find them with regex
        if len(matches) == 0:
            for regex_pattern in regex_patterns:
                matches += [*re.finditer(regex_pattern, allele_substring)]

            # We sort the matches, to replace the longest matching ones first.
            matches.sort(key=lambda match: len(match.group()), reverse=True)

        for match in matches:
            if match.group() in allele_substring:
                this_list = replace_substring_by_match(allele_substring, match)
                # Recursion
                this_list = replace_allele_features(
                    regex_patterns, this_list, matches)
                break
        else:
            # If none of the matches is in the allele_substring, we just return it as is.
            this_list = [allele_substring]

        out_list += this_list

    return out_list


def build_regex2syntax_rule(syntax_rules):
    """
    A dictionary in which the keys are the regex patterns and the values are the syntax_rules passed
    as arguments.
    """
    out_dict = dict()
    for syntax_rule in syntax_rules:
        out_dict[syntax_rule['regex']] = syntax_rule
    return out_dict


def sort_result(input_list: list[str, re.Match]) -> tuple[list[re.Match], list[str]]:
    """
    See docstring of `replace_allele_features` for relevant example. From the output of that function:

    input_list = [<re.Match for 'A'>, <re.Match for '321'>, <re.Match for 'B'>, '**']

    Returns two lists:
        * matches: list of re.matches only, in their appearing order
        * unmatched: list of unmatched strings, excluding those that are only made of symbols (matched by regex ^[^a-zA-Z\d]+$)

    For the example, it would return:
    matches = [<re.Match for 'A'>, <re.Match for '321'>, <re.Match for 'B'>]
    unmatched = [] # because '**' matches '^[^a-zA-Z\d]+$'
    """
    matches = list()
    unmatched = list()
    for r in input_list:
        # This is a match
        if type(r) != str:
            matches.append(r)
        # This is a string that has letters or digits in it, so we can't skip it
        elif not re.match('^[^a-zA-Z\d]+$', r):
            unmatched.append(r)
    return matches, unmatched


def allele_is_invalid(allele_description, syntax_rules, allele_type, allowed_types, gene):

    regex2syntax_rule = build_regex2syntax_rule(syntax_rules)

    result = replace_allele_features(
        list(regex2syntax_rule.keys()), [allele_description], [])

    # Filter out the non-digit non-letter characters
    matches, unmatched = sort_result(result)

    if len(unmatched):
        return 'pattern not matched\t' + ','.join(unmatched) + '\t'

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
        applied_rules.append(
            f'{syntax_rule["type"]}:{syntax_rule["rule_name"]}')

    error_list = list()
    if len(invalid_list):
        error_list.append(','.join(invalid_list))
    if len(sequence_error_list):
        error_list.append(','.join(sequence_error_list))
    if len(error_list):
        return 'invalid\t' + '|'.join(error_list)

    encountered_types = frozenset(encountered_types)
    correct_type = allowed_types[encountered_types]

    corrections = [[], []]
    if correct_type != allele_type:
        corrections[0].append('wrong type')
        corrections[1].append(f'change type to: {correct_type}')

    expected = ','.join(expected_list)
    if expected != allele_description:
        corrections[0].append('typo:')
        corrections[1].append(f'fix to: {expected}')

    corrections.append(applied_rules)
    if len(corrections[0]):
        return '\t'.join([','.join(c) for c in corrections])

    return False
