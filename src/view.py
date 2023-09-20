INDEX = '''
{HEAD}
{CONFIG}
{ALIGNMENT}
{DISTANCE}
{PERMUTATION}
{TIME}
'''

HEAD = '''
 ********************************************************************
 * ARTEMIS (Version 20230828)                                       *
 * using ARTEM to Infer Sequence alignment                          *
 * Reference: TODO                                                  *
 * Please email comments and suggestions to dbohdan@iimcb.gov.pl    *
 ********************************************************************
'''

CONFIG = '''
Configuration:
r={r}
q={q}
rformat={rformat}
qformat={qformat}
rres={rres}
qres={qres}
rresneg={rresneg}
qresneg={qresneg}
rseed={rseed}
qseed={qseed}
saveto={saveto}
saveformat={saveformat}
saveres={saveres}
threads={threads}
matchrange={matchrange}
stepdiv={stepdiv}
nlargest={nlargest}
shift={shift}
permutation={permutation}
verbose={verbose}
'''

ALIGNMENT = '''
Name of structure r: {rName}:{rChain}
Name of structure q: {qName}:{qChain} (to be superimposed onto structure r)
Length of structure r: {rLength} residues
Length of structure q: {qLength} residues

Aligned length= {aliLength}, RMSD= {RMSD:6.2f}, Seq_ID=n_identical/n_aligned= {Seq_ID:4.3f}
TM-score= {rTMscore:6.5f} (normalized by length of structure r: L={rLength}, d0={r_d0:.2f})
TM-score= {qTMscore:6.5f} (normalized by length of structure q: L={qLength}, d0={q_d0:.2f})

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
{rAlignment}
{distances}
{qAlignment}
'''

PERMUTATION = '''
_____________________________________________________________________
Alignment with permutations:
Aligned length= {p_aliLength}
TM-score= {p_rTMscore:6.5f} (normalized by length of structure r)
TM-score= {p_qTMscore:6.5f} (normalized by length of structure q)
RMSD= {p_RMSD:6.2f}
Seq_ID=n_identical/n_aligned= {p_Seq_ID:4.3f}
'''

TIME = '''
#Total CPU time is {time_total:5.2f} seconds
'''
