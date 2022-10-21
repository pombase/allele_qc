# %%

old_alignment = "VAQCIKVTVIFLAQCVKVTVIFL"
new_alignment = "VAQCIKVT----AQCVKVTVIFL"

new_seq = new_alignment.replace('-', '')
old_seq = old_alignment.replace('-', '')

# old 151 -> new 162


def get_other_index(this_alignment, other_alignment, this_index):

    count_other = -1
    count_this = -1

    for i in range(len(this_alignment)):
        count_this += this_alignment[i] != '-'
        count_other += other_alignment[i] != '-'
        if count_this == this_index:
            if other_alignment[i] == '-':
                return None
            return count_other
    # The coordinate does not exist in old one (new one is longer)
    return None


for i in range(len(new_alignment)):
    old_i = get_other_index(new_alignment, old_alignment, i)
    print(i, old_i)
    if old_i is None:
        continue
    # print(new_seq[i], old_seq[old_i])
