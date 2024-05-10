def len_seq(seq):
    ''' str -> int
    Returns the length of the sequence
    '''
    return len(seq)

def GC_content(seq):
    ''' str -> float
    Assumption: seq is a nucleotide sequence
    Returns the GC content
    '''
    # GC : int ; the counter of G and C
    GC = 0
    # base : str ; a nucleotide
    for base in seq:
        if base == 'G' or base == 'C':
            GC+=1
    return GC/len(seq)

def rev_compl(seq):
    ''' str -> str
    Assumption: seq is a nucleotide sequence
    Returns the reverse complement sequence
    '''
    # nucls: Dict{str:str} ; dictionary of complement nucleotides
    nucls = {'A':'T','T':'A','G':'C','C':'G'}
    # rev_compl : str ; future reverse complement sequence
    rev_compl = ''
    # base : str ; bases of the sequence
    for base in seq:
        rev_compl = nucls[base] + rev_compl
    return rev_compl

def ORF_position(seq):
    ''' str -> List[int],List[int]
    Assumption: seq is a nucleotide sequence
    Returns the positions (starting from 0) of ATG and STOP (TGA, TAA, TAG)
    for a given strand
    '''
    # ATGs : List[int] ; a list of positions of ATG
    ATGs = []
    # STOPs : List[int] ; a list of positions of ATG
    STOPs = []
    # i : int ; index in a sequence
    for i in range(len(seq)-2):
        if seq[i:i+3]=='ATG':
            ATGs.append(i)
        elif seq[i:i+3] in {'TGA','TAA','TAG'}:
            STOPs.append(i)
    return ATGs, STOPs

def compare_seqs(seq1,seq2):
    ''' str*str -> int,int,int
    Assumption: seq1 and seq2 are aligned sequences
    Returns the number of matches, substitutions and indels
    '''
    # matches : int ; counter for matches
    # substs : int ; counter for substitutions
    # indels : int ; counter for indels
    matches, substs, indels = 0,0,0
    # i : int ; index in the sequence
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matches+=1
        # Note that gap '-' will never be aligned with another gap for 2 seq. align.
        elif seq1[i] == '-' or seq2[i] =='-':
            indels+=1
        else:
            substs+=1
    return (matches, substs, indels)

def edit_distance(seq1,seq2):
    ''' str*str -> int
    Assumption: seq1 and seq2 are aligned sequences
    Returns the edit distance (the sum of indels and substitutions)
    '''
    # compseqs : tuple(int) ; number of matches, substitutions and indels
    compseqs = compare_seqs(seq1,seq2)
    return compseqs[1]+compseqs[2]

def align_score(seq1,seq2, Wm, Ws, Wi):
    ''' str*str*float*float*float -> float
    Assumption: seq1 and seq2 are aligned sequences
    Returns the alignment score - the weighted sum of matches, substitutions and indels
    (Wm,Ws,Wi are corresponding weights)
    '''
    # compseqs : tuple(int,int,int) ; number of matches, substitutions and indels
    compseqs = compare_seqs(seq1,seq2)
    return Wm*compseqs[0]+Ws*compseqs[1]+Wi*compseqs[2]

def identity_ratio(seq1,seq2):
    ''' str*str -> int,int,int
    Assumption: seq1 and seq2 are aligned sequences
    Returns the identity ratio - the proportion of matches
    '''
    return compare_seqs(seq1,seq2)[0]/len(seq1)

def find_all_ORFs(seq, length_limit=48):
    ''' str*int -> dict[str:List[tuple(int)]]
    Assumption: seq is a nucleotide sequence, length_limit is a positive integer
    corresponding to a minimal of a potential ORF (48 by default)
    Returns a dictionary with directions '+' and '-' as keys,
    and list of tuples of the positions of the start and the end of
    corresponding ORFs
    Note: positions start from 0 and not from 1
    '''
    # ORFs : dict[str:List[tuple(int)]] ; dictionary with directions '+' and '-' as keys,
    # and list of tuples of the positions of the start and the end of corresponding ORFs
    ORFs = {'+':[],'-':[]}

    # dir : str ; direction of the sequence, + for forward and - for reverse
    for dir in '+','-':
        if dir == '-':
            # seq : str ; input sequence, reversed or not
            seq = rev_compl(seq)
        # starts : List[int] ; list of start positions
        starts =[]
        # starts : List[int] ; list of end positions
        ends = []
        # i : int ; a position in the sequence
        for i in range(len(seq)-2):
            # triplet : str ; read triplet of nucleotides
            triplet = seq[i:i+3]
            if triplet == 'ATG':
                starts.append(i)
            elif triplet in {'TGA','TAA','TAG'}:
                ends.append(i+2)
        # start : int ; a position of a start
        for start in starts:
            # end : int ; a position of a potential end
            for end in ends:
                # length : int ; length of the potential ORF
                length = end-start+1
                if length>0 and length%3==0:
                    if length>length_limit:
                        ORFs[dir].append((start,end))
                    break
    return ORFs

def read_genetic_code(path_to_gencode):
    ''' str -> Dict[str:str]
    Assumption: path_to_gencode leads to a file where each represents "{Codon} {Aminoacid}"
    Returns a dictionary {Codon:Aminoacid}
    '''
    # gencode : Dict[str:str] ; dictionary for genetic code {Codon:Aminoacid}
    gencode = {}
    # streamr : stream ; input file stream
    with open(path_to_gencode) as streamr:
        # line : str ; a read line of the file
        for line in streamr:
            if len(line.split())==2:
                # codon, aa : str ; codon and aminoacid read from genetic code file
                codon, aa = line.split()
                gencode[codon]=aa
    return gencode 

def translate(ORF,gencode):
    ''' str * Dict[str:str] -> str
    Assumption: ORF is an open reading frame, gencode is a dictionary {Codon:Aminoacid}
    Returns a translated ORF according to given genetic code
    '''
    # protein : str ; a protein sequence
    protein = ''
    # i : int ; index in the ORF
    i = 0
    while i<len(ORF):
        # codon : str
        codon = ORF[i:i+3]
        protein+=gencode[codon]
        i+=3
    return protein

# TEST_SEQ : str ; test sequence 
TEST_SEQ = 'GCGCATGCTTACATAGGCCTAACAAACGGCTTATGACTAG'

# TEST_ALIGNED_SEQ : List[str] ; sequences for the question 2
TEST_ALIGNED_SEQ = ['AGTTAGAG--TAGGGCCAGCCAGCATAGCA-----GGTT',
         'AGTTCCAGTATAGGCCCA---AGCAAAGCAGTACCGGTT']

assert len_seq('ATGC')==4
assert GC_content('ATGC')==0.5
assert rev_compl('ATGC')=='GCAT'
assert ORF_position('GATGCATTAA') == ([1],[7])
assert compare_seqs('CA-CGTGCTGACCCAACC','CAGCGCGCTGG--CAGCC') == (12, 3, 3)
assert edit_distance('CA-CGTGCTGACCCAACC','CAGCGCGCTGG--CAGCC') == 6
assert align_score('CA-CGTGCTGACCCAACC','CAGCGCGCTGG--CAGCC',+2,-1,-2) == 15
assert identity_ratio('CA-CGTGCTGACCCAACC','CAGCGCGCTGG--CAGCC') == 12/18
assert find_all_ORFs('GCATGCGAGCTATATGCATGGCCTACACGCTAAAGATGCAGATGCTAAATACGGAAGCAGACCGTA\
GCGCTTGATACATCTCGTAAGCATA',48)=={'+': [(2, 73), (13, 66), (17, 73)], '-': [(1, 81), (12, 68)]}
assert translate('ATGAAATGA',{'ATG':'M','AAA':'K','TGA':'*'}) == 'MK*'
