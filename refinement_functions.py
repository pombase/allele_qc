import re
from models import SyntaxRule


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
        if type(allele_substring) == re.Match:
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


def build_regex2syntax_rule(syntax_rules: list[SyntaxRule]) -> dict[str, SyntaxRule]:
    """
    A dictionary in which the keys are the regex patterns and the values are the syntax_rules passed
    as arguments.
    """
    out_dict = dict()
    for syntax_rule in syntax_rules:
        out_dict[syntax_rule.regex] = syntax_rule
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


def get_allele_parts_from_result(result):
    """The result parts, excluding non-digit non-letter characters."""
    allele_parts = list()
    for r in result:
        if type(r) == str:
            if not re.match('^[^a-zA-Z\d]+$', r):
                allele_parts.append(r)
        else:
            allele_parts.append(r.group())

    return allele_parts


def check_allele_description(allele_description, syntax_rules, allele_type, allowed_types, gene):
    """
    Use replace_allele_features to identify patterns based on syntax rules, then validate
    the content of those patterns based on the grammar rules, and return output. See the
    example from test_data/allele_expected_results.tsv
    """
    regex2syntax_rule = build_regex2syntax_rule(syntax_rules)

    result = replace_allele_features(list(regex2syntax_rule.keys()), [allele_description], [])

    allele_parts = get_allele_parts_from_result(result)

    # Extract the matched and unmatched elements
    matches, unmatched = sort_result(result)

    output_dict = {
        'allele_parts': '',
        'needs_fixing': True,
        'change_description_to': '',
        'rules_applied': '',
        'pattern_error': '',
        'invalid_error': '',
        'sequence_error': '',
        'change_type_to': ''
    }

    if len(unmatched):
        output_dict['pattern_error'] = ','.join(unmatched)
        return output_dict

    # By default empty strings
    allele_part_types = ['' for m in matches]
    correct_name_list = ['' for m in matches]
    invalid_error_list = ['' for m in matches]
    sequence_error_list = ['' for m in matches]
    rules_applied = ['' for m in matches]

    for i, match in enumerate(matches):

        syntax_rule = regex2syntax_rule[match.re.pattern]
        allele_part_types[i] = syntax_rule.type

        rules_applied[i] = f'{syntax_rule.type}:{syntax_rule.rule_name}'

        invalid_error_list[i] = syntax_rule.check_invalid(match.groups())
        if invalid_error_list[i]:
            continue

        sequence_error_list[i] = syntax_rule.check_sequence(match.groups(), gene)
        correct_name_list[i] = syntax_rule.apply_syntax(match.groups())

    encountered_types = frozenset(allele_part_types)
    correct_type = allowed_types[encountered_types]

    if correct_type != allele_type:
        output_dict['change_type_to'] = correct_type

    correct_name = ','.join(correct_name_list)
    if correct_name != allele_description and all(correct_name_list):
        output_dict['change_description_to'] = correct_name

    output_dict['rules_applied'] = '|'.join(rules_applied) if any(rules_applied) else ''
    output_dict['invalid_error'] = '|'.join(invalid_error_list) if any(invalid_error_list) else ''
    output_dict['sequence_error'] = '|'.join(sequence_error_list) if any(sequence_error_list) else ''
    output_dict['allele_parts'] = '|'.join(allele_parts) if any(allele_parts) else ''

    must_be_empty = ['pattern_error', 'invalid_error', 'sequence_error', 'change_description_to', 'change_type_to']
    output_dict['needs_fixing'] = any(output_dict[key] for key in must_be_empty)

    return output_dict


def seq_error_change_description_to(allele_name, sequence_error):
    """
    Apply the proposed coordinate change in sequence_error:

    allele_name: 'K35A,P90T', sequence_error: 'K35A>K36A|' ==> returns K36A,P90T

    If no `>` in sequence_error or there are more than 1 ==> returns ''
    """
    new_allele_parts = list()
    for allele_part, error_message in zip(allele_name.split(','), sequence_error.split('|')):
        if '>' not in error_message:
            new_allele_parts.append(allele_part)
        elif error_message.count('>') == 1:
            old_val, new_val = error_message.split('>')
            new_allele_part = allele_part.replace(old_val, new_val)
            if new_allele_part == allele_part:
                return ''
            new_allele_parts.append(new_allele_part)
        else:
            return ''

    return ','.join(new_allele_parts)


def split_multiple_aa(value, regex):
    """Split into single variants: VLP-123-AAA => ['V123A', 'L124A', 'P125A']"""

    groups = re.match(regex, value).groups()
    return [f'{aa1}{int(groups[1])+i}{aa2}' for i, (aa1, aa2) in enumerate(zip(groups[0], groups[2]))]


def join_multiple_aa(values):
    """Opposite of split_multiple_aa"""
    sorted_values = sorted(values, key=lambda x: int(re.search(r'\d+', x).group()))
    return ''.join(v[0] for v in sorted_values) + '-' + sorted_values[0][1:-1] + '-' + ''.join(v[-1] for v in sorted_values)
