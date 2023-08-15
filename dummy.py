from ctd_support import ctd_check_sequence, ctd_convert_to_normal_variant

print(ctd_check_sequence('CTD-Y1F(r1-r12)', 'SPBC28F2.12'))
print(ctd_check_sequence('CTD-Y2F(r1-r12),C3A', 'SPBC28F2.12'))

print(ctd_convert_to_normal_variant('CTD-delta', 'SPBC28F2.12'))
print(ctd_convert_to_normal_variant('CTD-delta(r1-r1)', 'SPBC28F2.12'))
print(ctd_convert_to_normal_variant('CTD-delta(r1-r2)', 'SPBC28F2.12'))
print(ctd_convert_to_normal_variant('CTD-delta(r1-r3)', 'SPBC28F2.12'))
print(ctd_convert_to_normal_variant('CTD-delta(r1-r10-2)', 'SPBC28F2.12'))