import sys
sys.path.append('Project_library')

import re
import basic_seq_analysis as SeqAnalysis

def get_EMBL_summary_as_dict(path):
    ''' str -> dict[str:List[str]]
    Assumption: takes a path to EMBL database file
    Returns a dictionary of lists for each CDS UniProtID, each list consists of : ID of the sequence of origin,
    organism species, direction, beginning and end of CDS and function of the corresponding protein
    Note:!!! Dictionary with UniProtID doesn't seem to work since there're less CDSs selected when
    I used dictionary than when I used a list of lists (cf Output_Files)...
    '''
    # embl_sum : dict[str:List[str]] ; dictionary of lists {UniProtID : [Sequence ID, Species, Direction, Begin, End, Function]}
    embl_sum ={}
    # is_CDS : bool ; flag, the current line is in the CDS text block
    is_CDS = False
    # streamr : stream ; input file stream
    with open(path, 'r') as streamr:
        # line : str ; a read line of the file
        for line in streamr:
            # was_CDS : bool ; flag, the previous line is in the CDS text block
            was_CDS = is_CDS
            # the current line is still in the CDS block if the previous line was
            # in CDS and the current line starts with FT with FOUR spaces (instead of normal 3)
            is_CDS = (re.match('FT    ',line)!=None) and was_CDS
            # if the prevous line was the last line of CDS block, save the data
            if was_CDS and not is_CDS :
                embl_sum[uniprot]=[id, species, dir, begin, end, function]
            if re.match('OS   ', line):
                # scpecies : str ; species name
                species = line.split()[1] + ' '+ line.split()[2]

            elif re.match('ID   ', line):
                # id : str ; ID of the sequence
                id = line.split()[1][:-1]

            # if CDS block lines are read :
            elif is_CDS:
                # since 'product' text block sometimes extends beyond one line :
                # was_product : bool; flag, the previous line was in 'product' block
                # is_product : bool ; flag, the current line is in 'product' block
                was_product = is_product
                if re.match('FT                   /product="', line):
                    # function : str; function of a corresponding protein
                    function = re.split('"|\n',line)[1]
                    is_product = True
                if is_product and not re.match('FT                   /',line):
                    function = function +' '+ re.split(' +|\n',line)[1]
                if was_product and re.match('FT                   /',line):
                    is_product = False

                if re.match('FT                   /db_xref="UniProtKB', line) :
                    # uniprot : str ; UniProt ID of a corresponding protein
                    uniprot = line.split(':')[1][:-2]

            elif re.match('FT   CDS', line) and not re.search('[<>]',line):
                is_CDS = True
                is_product = False
                if re.search('complement', line):
                    # dir : str ; the direction of CDS
                    dir = '-'
                else:
                    dir = '+'
                # begin, end : str ; beginning and end of CDS
                begin, end = re.search(r'[0-9]+\.\.[0-9]+', line).group().split('..')
    return embl_sum


def write_EMBL_summary_file_from_dict(embl_sum, output_file_name):
    ''' dict[str:List[str]] * str -> None
    Assumption: embl_sum is the summary of EMBL data base file of the format returned
    by get_EMBL_summary_as_dict() function (dictionary of lists for UniProtID of each CDS consisting of :
    ID of the sequence of origin,organism species, direction, beginning and end of CDS and function
    of the corresponding protein)
    Writes a text file with an EMBL summary formatted as a table
    '''
    # header : str
    header = '{:^12s} {:^24s} {:^3s} {:>8s} - {:8s} {:^10s} {:60s}'.format('EMBL ID','Species','dir', 'begin','end', 'Uniprot','Function')
    # line_template : str ; template for each line of the table
    line_template = '{:12s} {:24s} {:^3s} {:>8s} - {:8s} {:10s} {:60s}'
    # dashes : str
    dashes = '{:12s} {:24s} {:^3s} {:>8s}---{:8s} {:10s} {:60s}'.format(12*'-',24*'-',3*'-',8*'-',8*'-',10*'-',60*'-')
    # streamw : stream ; output file stream in 'writing' mode
    with open(output_file_name, 'w') as streamw:
        streamw.write(dashes+'\n'+header+'\n'+dashes+'\n')
        # uniprot : str ; UniProtID, a key of input dictionary
        for uniprot in embl_sum:
            # line : List[str] ; a 'line' of EMBL summary
            line = embl_sum[uniprot][:5]+[uniprot]+embl_sum[uniprot][5:]
            streamw.write(line_template.format(*line)+'\n')
        streamw.write(dashes)

def get_and_select_CDS_positions(path_to_embl, path_to_regex):
    ''' str * str -> dict[str:List[str]]
    Assumptions: takes a path to EMBL database file and a path to a file with
    regular expressions needed for selection of CDSs of interest : regex for species name,
    gene name and protein function ('product') on 1st, 2nd and 3rd lines respectively
    with a new line character after each line
    Returns a dictionary of lists for each of selected CDS UniProtID, each list consists of :
    ID of the sequence of origin, direction, beginning and end of CDS
    '''
    # CDS_positions : dict[str:List[str]] ; dictionary of lists {UniProtID : [Sequence ID, Direction, Begin, End]}
    CDS_positions = {}
    # is_CDS : bool ; flag, the current line is in the CDS text block
    # species_OK : bool ; flag, the block corresponds to the right species name
    # gene_OK : bool ; flag, the current CDS text block has a right gene name
    # product_OK : bool ; flag, the current CDS text block has a right function (aka product)
    is_CDS, species_OK, gene_OK, product_OK = False,False,False,False
    # streamr_regex : stream ; input regex file stream
    with open(path_to_regex, 'r') as streamr_regex:
        # species_names, gene_names, product : str ; regex for species name, gene name
        # and protein function (aka product)
        species_names, gene_names, product = streamr_regex.read().splitlines()
    # streamr : stream ; input file stream
    with open(path_to_embl, 'r') as streamr:
        # line : str ; a read line
        for line in streamr:
            # was_CDS : bool ; flag, the previous line is in the CDS text block
            was_CDS = is_CDS
            is_CDS = (re.match('FT    ',line)!=None) and was_CDS
            if was_CDS and not is_CDS and (species_OK and (gene_OK or product_OK)):
                CDS_positions[uniprot] = [id, dir, begin, end]
            if re.match('OS   '+species_names, line):
                species_OK = True
            elif re.match('OS   ', line):
                species_OK = False

            if re.match('ID   ', line):
                # id : str ; ID of a sequence
                id = line.split()[1][:-1]

            if species_OK and is_CDS :
                if re.match('FT                   /gene="'+gene_names, line) :
                    gene_OK = True
                elif re.match('FT                   /product="'+product, line):
                    product_OK = True
                elif re.match('FT                   /db_xref="UniProtKB.+"', line) :
                    # uniprot : str ; UniProt ID of a corresponding protein
                    uniprot = line.split(':')[1][:-2]

            if re.match('FT   CDS', line) and not re.search('[<>]',line) and species_OK:
                is_CDS = True
                gene_OK = False
                product_OK = False
                if re.search('complement', line):
                    # dir : str ; the direction of CDS
                    dir = '-'
                else:
                    dir = '+'
                # begin, end : str ; beginning and end of CDS
                begin, end = re.search(r'[0-9]+\.\.[0-9]+', line).group().split('..')

    return CDS_positions

def read_fasta(path):
    ''' str -> dict[str:str]
    Assumption: input file is of a Fasta format, ';' is considered as commentary
    Returns a dictionary {sequence ID : sequence}
    '''
    # fasta : dict[str:str] ; the dictionary {Sequence ID : Sequence}
    fasta = {}
    # seq : str ; sequence
    seq = ''
    # id : str ; sequence ID
    id = None
    # streamr : stream ; input file stream
    with open(path,'r') as streamr:
        # line : str ; a read line
        for line in streamr:
            if line.strip()[0] != ';':
                if re.match('>',line):
                    if id != None:
                        fasta[id] = seq
                    id = line[1:-1]
                    seq = ''
                elif line[-1]=='\n':
                    seq = seq + line[:-1]
                else:
                    seq = seq + line
    fasta[id] = seq
    return fasta

def extract_CDSs_seq_from_fasta(CDS_positions,fasta):
    ''' dict[str:List[str]] * dict[str:str] -> dict[str:str]
    Assumption: CDS_positions is an output of my get_and_select_CDS_positions() function
    (dictionary of lists for each of selected CDS UniProtID, where each list consists of :
    ID of the sequence of origin, direction, beginning and end of CDS) ;
    fasta is an output of my read_fasta() function - a dictionary {Sequence ID : Sequence}
    Returns a dictionary {UniProtID : Sequence}
    '''
    # extracted_seqs : dict[str:str]; a dictionary {UniProtID : Sequence}
    extracted_seqs = {}
    # uniprot : str ; UniProt ID, a key of CDS_positions dictionary
    for uniprot in CDS_positions:
        # id, dir, begin, end: str; ID of the sequence of origin,direction, beginning and end
        id, dir, begin, end = CDS_positions[uniprot]
        if dir == '+':
            extracted_seqs[uniprot] = fasta[id][int(begin)-1:int(end)]
        else:
            extracted_seqs[uniprot] = SeqAnalysis.rev_compl(fasta[id][int(begin)-1:int(end)])
    return extracted_seqs


#################################################################
##################### IF I USE LIST OF LISTS ####################
#################################################################
def get_EMBL_summary_as_listoflists(path):
    ''' str -> List[List[str]]
    Assumption: takes a path to EMBL database file
    Returns a list of lists for each CDS ,each list consists of : ID of the sequence of origin,
    organism species, direction, beginning and end of CDS, UniProtID and function of the corresponding protein
    '''
    # embl_sum : List[List[str]] ; list of lists [Sequence ID, Species, Direction, Begin, End,UniProtID, Function] for each CDS
    embl_sum =[]
    # is_CDS : bool ; flag, the current line is in the CDS text block
    is_CDS = False
    # streamr : stream ; input file stream
    with open(path, 'r') as streamr:
        # line : str ; a read line of the file
        for line in streamr:
            # was_CDS : bool ; flag, the previous line is in the CDS text block
            was_CDS = is_CDS
            # the current line is still in the CDS block if the previous line was
            # in CDS and the current line starts with FT with FOUR spaces (instead of normal 3)
            is_CDS = (re.match('FT    ',line)!=None) and was_CDS
            # if the previous line was the last line of CDS block, save the data
            if was_CDS and not is_CDS :
                embl_sum.append([id, species, dir, begin, end,uniprot, function])
            if re.match('OS   ', line):
                # species : str ; species name
                species = line.split()[1] + ' '+ line.split()[2]

            elif re.match('ID   ', line):
                # id : str ; ID of the sequence
                id = line.split()[1][:-1]

            # if CDS block lines are read :
            elif is_CDS:
                # since 'product' text block sometimes extends beyond one line :
                # was_product : bool; flag, the previous line was in 'product' block
                # is_product : bool ; flag, the current line is in 'product' block
                was_product = is_product
                if re.match('FT                   /product="', line):
                    # function : str; function of a corresponding protein
                    function = re.split('"|\n',line)[1]
                    is_product = True
                if is_product and not re.match('FT                   /',line):
                    function = function +' '+ re.split(' +|\n',line)[1]
                if was_product and re.match('FT                   /',line):
                    is_product = False

                if re.match('FT                   /db_xref="UniProtKB', line) :
                    # uniprot : str ; UniProt ID of a corresponding protein
                    uniprot = line.split(':')[1][:-2]

            elif re.match('FT   CDS', line) and not re.search('[<>]',line):
                is_CDS = True
                is_product = False
                if re.search('complement', line):
                    # dir : str ; the direction of CDS
                    dir = '-'
                else:
                    dir = '+'
                # begin, end : str ; beginning and end of CDS
                begin, end = re.search(r'[0-9]+\.\.[0-9]+', line).group().split('..')
    return embl_sum

def write_EMBL_summary_file_from_listoflists(embl_sum, output_file_name):
    ''' List[List[str]] * str -> None
    Assumption: embl_sum is the summary of EMBL data base file of the format returned
    by get_EMBL_summary_as_listoflists() function (a list of lists for each CDS ,each
    list consists of : ID of the sequence of origin, organism species, direction, beginning
    and end of CDS, UniProtID and function of the corresponding protein)
    Writes a text file with an EMBL summary formatted as a table
    '''
    # header : str
    header = '{:^12s} {:^24s} {:^3s} {:>8s} - {:8s} {:^10s} {:60s}'.format('EMBL ID','Species','dir', 'begin','end', 'Uniprot','Function')
    # line_template : str ; template for each line of the table
    line_template = '{:12s} {:24s} {:^3s} {:>8s} - {:8s} {:10s} {:60s}'
    # dashes : str
    dashes = '{:12s} {:24s} {:^3s} {:>8s}---{:8s} {:10s} {:60s}'.format(12*'-',24*'-',3*'-',8*'-',8*'-',10*'-',60*'-')
    # streamw : stream ; output file stream in 'writing' mode
    with open(output_file_name, 'w') as streamw:
        streamw.write(dashes+'\n'+header+'\n'+dashes+'\n')
        # line : List[str] ; a 'line' of EMBL summary
        for line in embl_sum:
            streamw.write(line_template.format(*line)+'\n')
        streamw.write(dashes)
