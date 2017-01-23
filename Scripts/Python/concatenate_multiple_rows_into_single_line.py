#TRINITY_GG_428_c0_g1_i1_orf1 PF13499.1 EF_hand_5
#TRINITY_GG_428_c0_g1_i1_orf1 PF00036.27 efhand
#TRINITY_GG_428_c0_g1_i1_orf1 PF13405.1 EF_hand_4
#TRINITY_GG_428_c0_g1_i1_orf1 PF13833.1 EF_hand_6
#TRINITY_GG_428_c0_g1_i1_orf1 PF13202.1 EF_hand_3
#TRINITY_GG_429_c0_g1_i1_orf1 PF00156.22 Pribosyltran
#TRINITY_GG_431_c5_g1_i1_orf1 PF00475.13 IGPD
#TRINITY_GG_461_c0_g1_i1_orf1 PF01208.12 URO-D
#TRINITY_GG_461_c0_g1_i1_orf1 PF12876.2 Cellulase-like

import re

out_lines = []
with open('file.txt', 'r') as f:
    key = None
    key_lines = []
    for line in f:
        m = re.match(r'^(\S+)\s(.+)$', line)
        k, v = m.group(1), m.group(2)
        if k != key:
            if key:
                out_lines.append('{0} {1}'.format(key, ' | '.join(key_lines)))
            key = k
            key_lines = [v]
        else:
            key_lines.append(v)
    else:
        if key:
            out_lines.append('{0} {1}'.format(key, ' | '.join(key_lines)))

with open('out.txt', 'w') as f:
    f.write('\n'.join(out_lines) + "\n")
