# Aldol figure 1 (aldol scheme), the script as first written, unchanged apart from review/SVG calls
from bondline import *
st = Style(); S = Scene(st)

def acetaldehyde():
    m = Mol(); c1 = m.atom("C", (0, 0)); m.attach(c1, 210, "C"); m.attach(c1, -30, "O", 2)   # H3C-CH=O
    # zigzag: CH3 lower-left... use simple 2-atom drawing
    return m
def acet():
    m = Mol(); me = m.atom("C", (0, 0)); c = m.attach(me, 30, "C"); m.attach(c, -30, "O", 2); return m

a1, a2 = acet(), acet()
# product: 3-hydroxybutanal  CH3-CH(OH)-CH2-CHO
p = Mol(); c4 = p.atom("C", (0, 0)); c3 = p.attach(c4, 30, "C"); p.attach(c3, 90, "OH")
c2 = p.attach(c3, -30, "C"); c1 = p.attach(c2, 30, "C"); p.attach(c1, -30, "O", 2)

S.caption(a1, "ethanal"); S.caption(a2, "ethanal"); S.caption(p, "3-hydroxybutanal")
S.row([a1, "+", a2, S.arrow(above="NaOH (aq)", below="5 °C"), p])
print(S.save("aldol.png")); S.check()
