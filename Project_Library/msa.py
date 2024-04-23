# Alignement par paires et alignement multiple de séquences biologiques
# Dérivé de https://github.com/zahrasalarian/Bioinformatics-Projects

from math import factorial
from copy import deepcopy

ALPHABET = 'ACDEFGHIKLMNPQRSTVWY*' # Valide pour ADN et protéines

# Matrice par défaut : match = +3, mismatch = -1
DEFT_MTX = { (aa1, aa2) : 3 if aa1 == aa2 and aa1 != '*'
                    else -1 if aa1 != '*' and aa2 != '*'
                    else  0 if aa1 == '*' and aa2 == '*'
                    else -4
                    for aa1 in ALPHABET for aa2 in ALPHABET }

DEFT_GOP = -10.5 # Ouverture de gap
DEFT_GEP = -6.5  # Extension de gap

def pair_align(seq0, seq1, sub=DEFT_MTX, gop=DEFT_GOP, gep=DEFT_GEP,
                                                               end_free=False):
    ''' str * str * float * float * Dict[str, str: float] -> str * str * float
        Hypothèses : seq0 et seq1 non vides, tous les tuples de caractères de
          seq0 et seq1, et '-', sont des clés de sub; gex < gop.
        Retourne, pour la matrice de scores sub, le coût d'ouverture de gap
          gop et le coût d'extension de gap gex, l'alignement optimal de seq0
          et seq1, ainsi que le score d'alignement. Si end_free est True,
          réalise un alignement "end-gap free".
    '''
    # gap_indices : List[int] ; position des indels sur seq1 alignée avec seq0
    global gap_indices

    # l0, l1 : int ; longueurs de seq0 et seq1
    l0 = len(seq0) # 2e indice : j, horizontal)
    l1 = len(seq1) # 1er indice : i, vertical

    # i, j : int ; indices pour, respectivement, seq1 et seq0

    ''' Création des matrice S, D, H, V, chacune remplie de 0 '''
    # S : List[List[int]] ; scores d'alignement (l1+1 lignes de l0+1 colonnes)
    # D : List[List[int]] ; scores en diagonale (l1+1 lignes de l0+1 colonnes)
    # H : List[List[int]] ; scores si gaps hor. (l1+1 lignes de l0+1 colonnes)
    # V : List[List[int]] ; scores si gaps ver. (l1+1 lignes de l0+1 colonnes)
    S, D, H, V = [], [], [], []
    for i in range(l1 + 1):
        S.append([0] * (l0 + 1))
        D.append([0] * (l0 + 1))
        H.append([0] * (l0 + 1))
        V.append([0] * (l0 + 1))

    ''' Remplissage de la première ligne et de la première colonne de S '''
    if not end_free:
        #S[0][0] = 0
        for i in range(1, l1 + 1):
            S[i][0] = S[i - 1][0] + gep
            H[i][0] = gop
            D[i][0] = gop
        for j in range(1, l0 + 1):
            S[0][j] = S[0][j - 1] + gep
            V[0][j] = gop
            D[0][j] = gop

    ''' Création de la matrice O, remplie de "☺" '''
    # O : List[List[char]] ; orientations (l1+1 lignes de l0+1 colonnes)
    O = []
    for i in range(l1 + 1):
        O.append(['☺'] * (l0 + 1))

    ''' Remplissage de la première ligne et de la première colonne de O '''
    for i in range(1, l1 + 1):
        O[i][0] = '↑'
    for j in range(1, l0 + 1):
        O[0][j] = '←'

    ''' Remplissage des matrices de score et de direction '''
    for i in range(1, l1 + 1):
        for j in range(1, l0 + 1):
            if seq1[i - 1] == '-' and seq0[j - 1] == '-':
                D[i][j] = S[i - 1][j - 1]
            elif seq1[i - 1] == '-' or seq0[j - 1] == '-':
                D[i][j] = S[i - 1][j - 1] + gep
            else:
                # score pour ajouter un (més)appariement
                D[i][j] = S[i - 1][j - 1] + sub[seq1[i - 1], seq0[j - 1]]
            # score pour ajouter un gap horizontal
            H[i][j] = max(D[i][j - 1] + gop, H[i][j - 1] + gep)
            # score pour ajouter un gap vertical
            V[i][j] = max(D[i - 1][j] + gop, V[i - 1][j] + gep)
            # smax : float ; score optimal jusqu'ici
            smax = D[i][j]
           # odir : str ; direction du chemin optimal en i, j
            odir = '↖'
            if H[i][j] > smax:
                smax = H[i][j]
                odir = '←'
            if V[i][j] > smax:
                smax = V[i][j]
                odir = '↑'
            S[i][j] = smax
            O[i][j] = odir

    ''' Affichage des matrices '''
    if False: # if True pour afficher ; if False pour ne pas afficher
        for i in range(l1+1):
            for j in range(l0+1):
                print('  {:s} {:+03d}'.format(O[i][j],S[i][j]),end='')
            print()

    ''' Création de l'alignement '''
    gap_indices = []
    # align0, align1 : str ; respectivement seq0 et seq1 alignées
    align0 = ''
    align1 = ''
    j = l0
    i = l1

    while i > 0 or j > 0:
        if O[i][j] == '↖':
            if seq0[j - 1] != '-' or seq1[i - 1] != '-':
                align0 = seq0[j - 1] + align0
                align1 = seq1[i - 1] + align1
            i = i - 1
            j = j - 1
        elif O[i][j] == '←' :
            align0 = seq0[j - 1] + align0
            align1 = '-' + align1
            if len(gap_indices) != 0:
                gap_indices = [x+1 for x in gap_indices]
                gap_indices.append(i)
            else:
                gap_indices.append(j)
            j = j - 1
        else:
            align0 = '-' + align0
            align1 = seq1[i - 1] + align1
            i = i - 1

    if False:
        print('\n',seq0,'\n',seq1,'\n',gop,gep,'->',S[l1][l0],'\n',
              align0,'\n',align1,'\n')
    return align0, align1, S[l1][l0]

def calc_MSA_score(seqs, sub=DEFT_MTX, gop=DEFT_GOP, gep=DEFT_GEP):
    ''' List[str] * Dict[str,str:float] * float * float -> float
        Hypothèse : toutes les chaînes de seqs ont la même longueur.
        Retourne le score de l'alignement multiple seqs, étant donnés les
          scores dans la matrice sub, et les poids gop et gep pour les indels.
    '''
    # score : int ; le score, à retourner
    score = 0
    start_gap = False
    # i : int ; indice sur la position dans l'alignement
    for seq in seqs.values():
        lseq = len(seq)
        break
    for i in range(lseq):
        # j, k : int ; indices de deux séquences
        for seqj in seqs:
            for seqk in seqs:
                if seqj <= seqk:
                    continue
                if (seqs[seqj][i], seqs[seqk][i]) in sub:
                    score += sub[seqs[seqj][i], seqs[seqk][i]]
                elif (seqs[seqj][i]=='-' and (i==0 or seqs[seqj][i-1]!='-')) or\
                     (seqs[seqk][i]=='-' and (i==0 or seqs[seqk][i-1]!='-')):
                    score += gop
                else:
                    score += gep
    return score

def calc_MSA_seqs(s2align, sub=DEFT_MTX, gop=DEFT_GOP, gep=DEFT_GEP):
    ''' List[str] * Dict[str,str:float] * float * float -> float
        Calcule l'alignement multiple optimal pour les séquences de s2align,
          en utilisant la matrice sub pour les substitutions, et les poids gop
          et gep, respectivement pour les ouvertures et extensions d'indels.
    '''
    global gap_indices

    # Crée la matrice de scores par paires
    pairwise_similarities_matrix = {}
    for seqi in s2align:
        for seqj in s2align:
            if seqi == seqj:
                pairwise_similarities_matrix[(seqi,seqj)] = 0
            else:
                pairwise_similarities_matrix[(seqi,seqj)] = \
                     pair_align(s2align[seqi], s2align[seqj], sub, gop, gep)[2]

    # Additionne les scores pour chaque séquence
    sum_score_seqs = {}
    for seqi in s2align:
        sum_score_seqs[seqi] = \
           sum([pairwise_similarities_matrix[(seqi,seqj)] for seqj in s2align])

    # Trouve la séquence centrale = celle dont la somme des scores est maximum
    center_key = max(sum_score_seqs, key=sum_score_seqs.get)

    # Trie les séquences par ressemblace avec la centrale
    seqs_in_descending_order = \
        sorted([(seqj,pairwise_similarities_matrix[(center_key,seqj)])
                        for seqj in s2align], key=lambda k: k[1], reverse=True)
    seqs_in_descending_order.remove((center_key,0))
    #print(center_key, seqs_in_descending_order)

    # Aligne chaque séquence à la séquence centrale
    alignments_with_center = {}
    center_seq = s2align[center_key]
    for (seqi,_) in seqs_in_descending_order:
        align_X, align_Y, s = pair_align(s2align[seqi], center_seq, sub,
                                                                      gop, gep)
        alignments_with_center[seqi] = [align_X, align_Y, s]
        center_seq = align_Y #center_seq = deepcopy(align_Y)
        for (seqj,_) in seqs_in_descending_order:
            if seqj == seqi:
                break
            if seqj != center_key:
                gap_indices.reverse()
                for k in gap_indices:
                    alignments_with_center[seqj][0] = \
                    alignments_with_center[seqj][0][:k] + '-' + \
                    alignments_with_center[seqj][0][k:]
                alignments_with_center[seqj] = \
                    [alignments_with_center[seqj][0], center_seq,
                     alignments_with_center[seqj][2]]
    MSA_seqs = {}
    MSA_seqs[center_key] = center_seq
    for k, v in alignments_with_center.items():
        MSA_seqs[k] = v[0]

    return MSA_seqs, calc_MSA_score(MSA_seqs, sub, gop, gep)

def pick_blocks(seqs):
    ''' List[str] -> List[i]
        Récupère les positions des blocs à partir de seqs.
    '''
    seq_num = len(seqs)
    blocks = []
    for seq in seqs.values():
        lseq = len(seq)
        break
    for i in range(lseq):
        matches_num = 0
        for seqj in seqs:
            for seqk in seqs:
                if seqj <= seqk:
                    continue
                if seqs[seqj][i] == seqs[seqk][i]:
                    matches_num += 1
        if matches_num != factorial(seq_num)/((2)*factorial(seq_num - 2)):
            if len(blocks) % 2 == 0:
                blocks.append(i)
        elif matches_num == factorial(seq_num)/((2)*factorial(seq_num - 2)):
            if len(blocks) > 0:
                if i - blocks[-1] > 1:
                    blocks.append(i)
                else:
                    blocks.pop()
    if len(blocks) % 2 == 1:
        blocks.append(i)
    return blocks

def delete_gaps_from_block(MSA_seqs, start, end):
    ''' List[str] * int * int -> List[str]
        Supprime les gaps de toutes les séquences, entre start et end
    '''
    MSA_seqs = deepcopy(MSA_seqs)
    #print(MSA_seqs)
    for seqi in MSA_seqs:
        MSA_seqs[seqi] = MSA_seqs[seqi][start:end+1].replace('-', '')
    #print(MSA_seqs)
    return MSA_seqs

def star_align(seqs=None, sub=DEFT_MTX, gop=DEFT_GOP, gep=DEFT_GEP):
    ''' List[str] * Dict[str,str:float] * float * float -> float * List[str]
        Réalise l'alignement multiple optimal en utilisant l'algorithme
        d'alignement en étoile.
    '''
    if seqs is None:
        print('Erreur : pas de séquences à aligner')
        exit(-1)
    seq_num = len(seqs)
    gap_indices = []
    MSA_seqs, MSA_score = calc_MSA_seqs(seqs, sub, gop, gep)
    MSA_score_2 = None
    while MSA_score_2 is not None and MSA_score_2 > MSA_score:
        blocks = pick_blocks(MSA_seqs)
        for i in range(0, len(blocks), 2):
            start, end = blocks[i], blocks[i+1]
            if start == end:
                continue
            block_without_gap = delete_gaps_from_block(MSA_seqs, start, end)
            MSA_seqs_2, MSA_score_2 = calc_MSA_seqs(block_without_gap, sub,
                                                                      gop, gep)
            for i in range(seq_num):
                MSA_seqs_2[i] = MSA_seqs[i][:start] + MSA_seqs_2[i] +\
                                                            MSA_seqs[i][end+1:]
            MSA_score_2 = calc_MSA_score(MSA_seqs_2, sub, gop, gep)
            if MSA_score_2 > MSA_score:
                MSA_seqs = deepcopy(MSA_seqs_2)
                MSA_score = MSA_score_2
    seqout = {}
    # Remet les séquences dans l'ordre initial /!\ TODO : choix ini ou ali
    for seq in list(seqs.keys()):
        seqout[seq] = MSA_seqs[seq]
    return MSA_score, seqout

if __name__=="__main__":

    # Séquences de test
    DEFT_SEQ = {'SeqA':'MSAGAEKQDRWIALHG',
                'SeqB':'MSNGAEQERWIALVDHG',
                'SeqC':'MSNAKGERLALVCHG',
                'SeqD':'MSAVHKQDRLALVDIG'}
    MSA_score, MSA_seqs = star_align(DEFT_SEQ, DEFT_MTX, DEFT_GOP, DEFT_GEP)
    print(f'{len(DEFT_SEQ)} séquences, score = {MSA_score}')
    for seq in MSA_seqs:
        print(seq, MSA_seqs[seq])
