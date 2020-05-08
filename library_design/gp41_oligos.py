"""Script for making oligos for gp41 'mutational scan' phage display assay."""

import Bio.Seq
import Bio.SeqUtils.CodonUsageIndices
import Bio.Alphabet.IUPAC
import itertools
import pandas
import os

from Bio.Seq import Seq
from operator import itemgetter

viruses = ['BG505', 'BF520','ZA1197'] #List of file names for creating libraries
oligo_length = 93
tile = 3

avoidmotifs = ['GAATTC', 'AAGCTT']

linker5 = 'GGTGGTGGAGGTTCCGGGGGAGGAGGTTCGGGCGGTGGGGGAAGT'
linker3 = 'GGTGGTGGAGGTTCCGGGGGAGGAGGTTCGGGCGGTGGGGGAAGT'
adaptor5 = 'aggaattctacgctgagt' 
adaptor3 = 'tgatagcaagcttgcc' 

cai = Bio.SeqUtils.CodonUsageIndices.SharpEcoliIndex

translation_table = {}
codons = []
amino_acids = []
rt_table = {}

for (_nt1, _nt2, _nt3) in itertools.product(Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters, repeat = 3):
    codon = _nt1 + _nt2 + _nt3
    codons.append(codon)
for codon in codons:
    translation_table[codon] = str(Bio.Seq.Seq(codon).translate())
    if translation_table[codon] not in amino_acids:
        amino_acids.append(translation_table[codon])
        rt_table[translation_table[codon]] = [codon]
    else:
        rt_table[translation_table[codon]].append(codon)

score_dict = {} #Dict of all syn. codons & cai scores for an aa. No stop codons.
for amino_acid in amino_acids:
    if amino_acid != '*':
        scores = []
        for codon in rt_table[amino_acid]:
            scores.append((codon, cai[codon]))
        score_dict[amino_acid] = scores


def translate(oligo):
    """Script for translating DNA more quickly than Bio.Seq

    Args:
        `oligo` (str)
            DNA sequence to translate

    Returns:
        `aa_read` (str)
            Translated DNA sequence

    >>> oligo = 'ATGGGC'
    >>> translate(oligo) == 'MG'
    True
    """
    assert len(oligo) % 3 == 0, "Not translatable due to length"
    oligo = oligo.upper()
    aa_read = ''.join([translation_table[oligo[i : i + 3]] for i in 
            range(0, len(oligo), 3)])
    return aa_read


def get_seq(seq_file):
    """Return sequence from a `.fasta` file `seq_file`
    """
    seq_by_line = []
    with open(seq_file) as f:
        for line in f:
            if '>' in line:
                continue # Can input sequence with fasta header
            else:
                seq = line.strip()
                seq_by_line.append(seq)
    seq = ''.join(seq_by_line)
    return seq


def subseq_rand_mid(full_seq):
    """Subsequence viral seqs into oligos of *oligo_length* tiled by *tile*.

    The middle amino acid in each oligo is randomized so that each region of 
    gp41 is covered by 20 oligos each with a different amino acid at the 
    middle location.

    Args: 
        `full_seq` (str)
            Full DNA sequence (incl linkers) to be split up into oligos.

    Returns:
        `subseqs` (list)
            List of subsequenced oligos with randomized middle amino acid.
    """
    subseqs = []
    i = 0
    while i <= (len(full_seq) - oligo_length):
        subseq_front = full_seq[i : i + ((oligo_length-3)/2)]
        subseq_end = full_seq[i + (((oligo_length-3)/2)+3) : i + oligo_length]
        for aa in amino_acids:
            if aa != '*':
                subseq_total = subseq_front + max(score_dict[aa], 
                        key=lambda item:item[1])[0] + subseq_end
                assert len(subseq_total) == oligo_length, 'Oligos not correct length'
                subseqs.append(subseq_total)
        i = i + tile
    return subseqs
    

def remove_rsites(oligos, codon_scores_by_aa):
    """Remove restriction sites from subsequenced oligos.

    Args:
        `oligos` (list)
            List of oligos. Will remove rsites from each oligo.

        `codon_scores_by_aa` (dict)
            Dictionary of all codons and cai scores for each amino acid.
                e.g. {'C': [('TGT', 0.5), ('TGC', 1)]}
    
    Returns:
        `oligos` (list)
            Returns the list of oligos with restriciton sites removed. 
    """
    rsites_count = 0
    replace = False
    for motif in avoidmotifs:
        for n in range(len(oligos)):
            oligos[n] = oligos[n].upper()
            if motif in oligos[n]:
                aa_withrs = translate(oligos[n])
                replace = True
                while replace:
                    for i in range(len(oligos[n]) - len(motif) + 1):
                        seq = oligos[n][i : i + len(motif)] # sequence starting at i
                        if seq == motif:
                            done = False #need to replace a codon in this seq
                            rsites_count += 1
                            if oligos[n][i-(i%3):i-(i%3)+3] == 'ATG' or oligos[n][i-(i%3):i-(i%3)+3] == 'TGG': #Make sure codon to replace has synonymous codons
                                j = i + 3
                                print('Codon to replace ({0}) has no synonymous codons. Replace next codon: {1}'). \
                                        format(oligos[n][i-(i%3):i-(i%3)+3], oligos[n][j-(j%3):j-(j%3)+3])
                                i = j
                            for r in codon_scores_by_aa: #For each amino acid...
                                for codon_tup in codon_scores_by_aa[r]: #...And for each (codon, score) tuple that encodes that aa...
                                    if not done: #If codon at this location has not been replaced
                                        if oligos[n][i-(i%3):i-(i%3)+3] in codon_tup: #If codon we want to replace (first codon that contains part of restriciton site) is in score tuple codon_tup...
                                            l = codon_scores_by_aa[r] #...Turn all synonymous codons and their scores into a list
                                            l.sort(key=itemgetter(1), 
                                                    reverse=True)
                                            if oligos[n][i-(i%3):i-(i%3)+3] == l[0][0]:
                                                print('Replacing highest scoring codon ({0}) with '\
                                                        '2nd highest scoring: {1}').format(l[0][0], l[1][0])
                                                new_codon = l[1][0]
                                            elif oligos[n][i-(i%3):i-(i%3)+3] != l[0][0]:
                                                print('Replacing codon ({0}) with highest scoring codon: {1}').\
                                                        format(oligos[n][i-(i%3):i-(i%3)+3], l[0][0])
                                                new_codon = l[0][0]
                                            clean_oligo = oligos[n][:i-(i%3)] \
                                                    + new_codon \
                                                    + oligos[n][i-(i%3)+3:]
                                            oligos[n] = clean_oligo
                                            done = True #Codon has been replaced, don't look at other possible synonymous codons.
                    if motif not in oligos[n]: #Make sure motif not in oligo after going through replacements
                        replace = False
                assert len(oligos[n]) == oligo_length, 'Oligo lengths not maintained ' \
                        'while removing restriction sites.'          
                assert aa_withrs == translate(oligos[n]), 'The amino acid sequence ' \
                        'has been altered by removing restriction site.'
    print('{0} restriction sites were removed with synonymous substitution.'
            .format(rsites_count))
    return oligos


def main():

    virus_oligos = {}

    for virus in viruses:
        print('\nDesigning oligos for {0}.'.format(virus))
        opt_seq_file = 'sequences/{0}_gp41_dna_opt.fasta'.format(virus)
        opt_dna_seq = get_seq(opt_seq_file)

        assert translate(opt_dna_seq).upper() == get_seq('sequences/{0}_gp41_prot.fasta'
                .format(virus)).upper(), 'Translated seq does not equal given protein sequence.'

        full_seq = linker5 + opt_dna_seq + linker3
        assert len(full_seq) % 3 == 0, "Adding linkers made sequence untranslatable"
        
        subsequences = subseq_rand_mid(full_seq)

        cleaned_oligos = remove_rsites(subsequences, score_dict)

        oligos_with_translation = []
        oligo_count = 0
        for oligo in cleaned_oligos:
            translated_oligo = translate(oligo)
            random_aa = translated_oligo[(len(translated_oligo)-1)/2]
            final_oligo = adaptor5 + oligo + adaptor3
            oligos_with_translation.append((final_oligo, translated_oligo, random_aa, oligo_count))
            oligo_count += 1 #Keep track of oligo number separately for each virus
        
        virus_oligos[virus] = oligos_with_translation

    all_oligos_df = pandas.melt(pandas.DataFrame(virus_oligos), value_vars=viruses, 
            var_name='Virus', value_name='Oligo_Info')

    all_oligos_df[['Oligo', 'Prot', 'Rand_AA', 'Count']] = all_oligos_df['Oligo_Info'].apply(pandas.Series)
    print('{0} oligos designed.').format(len(all_oligos_df))
    all_oligos_df = all_oligos_df.drop(['Oligo_Info'], axis=1).drop_duplicates(['Prot'], keep='first')
    print('{0} oligos remain after dropping duplicate protein sequences.').format(len(all_oligos_df))

    all_oligos_df['Rand_Loc'] = all_oligos_df['Count'] // 20 + 1
    all_oligos_df = all_oligos_df.reset_index(drop=True)

    if not os.path.isdir('./oligos'):
        os.mkdir('./oligos')
        
    all_oligos_df.to_csv('./oligos/gp41_final_oligos.txt', index_label='Oligo_Num', 
            columns=['Virus', 'Rand_Loc', 'Rand_AA', 'Oligo'])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()