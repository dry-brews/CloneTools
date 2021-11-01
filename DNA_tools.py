import sys
import random

def reverse_complement(DNA_seq, type = "DNA"):
    if type == "DNA":
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
                      'a':'t', 't':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    elif type == "RNA":
        complement = {'A':'U', 'U':'A', 'G':'C', 'C':'G',
                      'a':'u', 'u':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    else:
        sys.stderr.write("Did not regognize sequence type: %s" % type)
        sys.exit()
    revcomp = []
    for i in range(len(DNA_seq))[::-1]:
        try:
            revcomp.append(complement[DNA_seq[i]])
        except KeyError:
            sys.stderr.write("Did not recognize base: %s" % DNA_seq[i])
            sys.exit()
    return ''.join(revcomp)

def align_primers(primer1, primer2, direction = "revcomp", min_overlap = 10):
    hit = False
    if direction == "revcomp":
        target = "-" * (len(primer1)-min_overlap) + reverse_complement(primer2) + "-" * (len(primer1)-min_overlap)
    for i in range(len(target)):
        hitscore = 0
        for j in range(len(primer1)):
            try:
                if primer1[j].upper() == target[i+j].upper():
                    hitscore +=1
                    continue
            except IndexError:
                break
            if target[i+j] == "-":
                continue
            else:
                hitscore = -100
                break
        if hitscore > min_overlap:
            hit = True
            query = "-" * i + primer1 + "-" * (len(target) - i - len(primer1))
            print_query = []
            print_target = []
            for i in range(len(query)):
                if query[i] != "-" or target[i] != "-":
                    print_query.append(query[i])
                    print_target.append(target[i])
            print(''.join(print_query))
            print(''.join(print_target))
    if hit == False:
        print("No suitable alignment found")
    return hit


def PCR(template, primer1, primer2, min_anneal_len = 14):
    template = template.upper()
    templateR = reverse_complement(template).upper()
    prime_ann1 = primer1.upper()[0-min_anneal_len:]
    prime_ann2 = primer2.upper()[0-min_anneal_len:]
    seeds = False
    if prime_ann1 in template and prime_ann2 in templateR:
        seeds = "fwd"
    elif prime_ann1 in templateR and prime_ann2 in template:
        seeds = "rev"
    if seeds == False:
        print(prime_ann1)
        print(prime_ann2)
        sys.stderr.write("couldn't find binding sites for primers\n")
        sys.exit()
    elif seeds == "fwd":
        primer2 = reverse_complement(primer2)
    elif seeds == "rev":
        temp = primer2
        primer2 = reverse_complement(primer1)
        primer1 = temp
    index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
    index_end = template.find(primer2.upper()[:min_anneal_len])
    if index_start < index_end:
        return(primer1 + template[index_start:index_end] + primer2)
    elif index_start > index_end:
        template = template[(index_start + index_end)//2:] + template[:(index_start + index_end)//2]
        index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
        index_end = template.find(primer2.upper()[:min_anneal_len])
        return(primer1 + template[index_start:index_end] + primer2)
    
def read_fasta(fasta_file):
    with open(fasta_file, "r") as file_in:
        read_next_line = False
        seq_list = []
        for line in file_in:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                read_next_line = True
            elif read_next_line == True:
                if line.startswith(">"):
                    sys.stderr.write("There should only be one sequence in the fasta. Exiting...")
                    sys.exit()
                seq_list.append(line.strip())
        file_in.close()
    seq = ''.join(seq_list)

    return (seq_name, seq)

def align_seqs(seq1, seq2, min_overlap = 0.8):
    #Quick alignment with no gaps for finding offsets between near-identical seqs
    #Returns offset to give more than min_overlap identities between seqs
    
    #query is the first seq
    query = seq1
    #set minimum aligned length as a fraction of the query length
    min_overlap = int(min_overlap * len(query))
    #pad target seq with '-' on both sides to slide the query along as a window
    target = "-" * (len(query)-min_overlap) + seq2 + "-" * (len(query)-min_overlap)
    hit = None

    for i in range(len(target) - len(query) - 1):
        hitscore = 0
        for j in range(len(query)):
            if target[i+j] == "-":
                continue
            if query[j] == "-":
                continue
            elif query[j].upper() == target[i+j].upper():
                hitscore +=1
                continue
        #if you find a place where >min_overlap residues match perfectly...
        if hitscore > min_overlap:
            #We are at the position i where query[j] == target[i+j]
            #If the seqs are already aligned, offset should be zero, but...
            #i will be (len(query)-min_overlap). Therefore offset = i - (len(query)-min_overlap)
            hit = i - (len(query)-min_overlap)
            break
    #if hit == "False":
    #    print("No suitable alignment found")
    return hit

def assemble_fragments_PCA(list_of_subseqs, min_anneal_len = 15, max_anneal_len = 30):
    #Takes linear DNA fragments with end homology to other sequences and assembles
    #them into a single piece of DNA. Works with short homology (i.e., GoldenGate)
    #or longer homology (i.e., PCA). Fragments need not be in order, but if one
    #fragment has homology to multiple others, it will likely cause issues

    frag_map = {} #Each fragment points to the the next fragment in the assembly
    offset_map = {} #Each fragment tells you how much to clip off the end before appending the next

    #look at every pair of fragments
    for f1 in list_of_subseqs:
        for f2 in list_of_subseqs: 
            #try to 3' end of f1 to the 5' end of f2 (fast gapless alignment)
            hit = align_seqs(
                f1[-max_anneal_len:],
                f2[:max_anneal_len],
                min_overlap = min_anneal_len / max_anneal_len
                )
            #Usually the alignment will fail and return None
            if hit != None:
                #If it does hit, point f1 to f2 and record the overlap / offset
                frag_map[f1] = f2
                offset_map[f1] = (hit - max_anneal_len)

    #Now, you have the alignments, go from one fragment to the next in line
    #and add each to the assembly (clipped based on the overlap / offset)

    f = sorted(frag_map.keys())[0]
    frag1_len = len(f)
    assembly = ''
    for i in range(len(frag_map)):
        try:
            assembly += f[:offset_map[f]]
        #the last frag in frag_map doesn't point to anything and doesn't have an offset_map entry
        except KeyError:
            assembly += f

    #maybe the assembly is circular. In this case, the first frag will appear at
    #both the beginning and end of the assembly. Try to align the 5' end of the
    #assembly to the 3' end of the assembly
    close_it = align_seqs(
        assembly[:frag1_len],
        assembly[-frag1_len:],
        min_overlap = 1. - (min_anneal_len / frag1_len)
        )

    #If it does not appear circular, return the linear assembly
    if close_it == None:
        return {
            "sequence": assembly,
            "topology": "linear",
            "strandedness": "ds",
            "sticky_ends": None
            }
    #If it does appear circular, clip off the repeat and return a circular piece of DNA
    else:
        return {
            "sequence": assembly[:frag1_len - close_it],
            "topology": "circular",
            "strandedness": "ds",
            "stick_ends": None}

RE_lib = {
    "SapI" : {
        "recog_site": "GCTCTTC",
        "cut1_pos": 8,
        "cut2_pos": 11,
        "sticky_ends": "5_prime_overhang"
    },
    
    "BsaI" : {
        "recog_site": "GGTCTC",
        "cut1_pos": 7,
        "cut2_pos": 11
    },
    
    "PaqCI" : {
        "recog_site": "CACCTGC",
        "cut1_pos": 11,
        "cut2_pos": 15
    },
    
    "BfuAI" : {
        "recog_site": "ACCTGC",
        "cut1_pos": 10,
        "cut2_pos": 14
    },
    }
def TypeIIsDigest(sequence, enzyme = "SapI"):
    recog_site=RE_lib[enzyme]["recog_site"]
    cut1_pos=RE_lib[enzyme]["cut1_pos"]
    cut2_pos=RE_lib[enzyme]["cut2_pos"]
    
    sequence = reverse_complement(sequence)
    if sequence.find(recog_site) != -1:
        cutsite1 = sequence.find(recog_site)+cut1_pos
        cutsite2 = sequence.find(recog_site)+cut2_pos
        sequence = sequence[cutsite1:cutsite2].lower() + sequence[cutsite2:]
    sequence = reverse_complement(sequence)
    if sequence.find(recog_site) != -1:
        cutsite1 = sequence.find(recog_site)+cut1_pos
        cutsite2 = sequence.find(recog_site)+cut2_pos
        sequence = sequence[cutsite1:cutsite2].lower() + sequence[cutsite2:]
 
    return sequence

def GoldenGateAssembly(list_of_seqs, enzyme = "SapI"):
    list_of_cut_seqs = [
        TypeIIsDigest(seq, enzyme=enzyme) for seq in list_of_seqs
    ]

    frag_map = {} #Each fragment points to the the next fragment in the assembly
    sticky_len = RE_lib[enzyme]["cut2_pos"] - RE_lib[enzyme]["cut1_pos"]

    #look at every pair of fragments
    for f1 in list_of_cut_seqs:
        for f2 in list_of_cut_seqs: 
            #try to 3' end of f1 to the 5' end of f2 (exact match only)
            if f1[-sticky_len:] == f2[:sticky_len]:
                #check that you've got lowercase on both ends
                if f1[-sticky_len:] == f1[-sticky_len:].lower() and f2[:sticky_len] == f2[:sticky_len].lower():
                    frag_map[f1] = f2

    #Now, you have the alignments, go from one fragment to the next in line
    #and add each to the assembly (clipped based on the sticky end len)

    f = sorted(frag_map.keys())[0]
    f_init = f
    assembly = f
    for i in range(len(frag_map)):
        # go to the next fragment
        f = frag_map[f]
        if f == f_init and i != 0:
            break
        #clip off the sticky end
        assembly = assembly[:-sticky_len]
        #and add the next fragment
        assembly += f
        #the last frag in frag_map doesn't point to anything and doesn't have an offset_map entry


    #maybe the assembly is circular. In this case, you'll have lowercase sticky
    #ends on both the start and the end
    if assembly[:sticky_len] == assembly[-sticky_len:]:
        if assembly[:sticky_len] == assembly[:sticky_len].lower() and assembly[-sticky_len:] == assembly[-sticky_len:].lower():
            return {
            "sequence": assembly[:-sticky_len].upper(),
            "topology": "circular",
            "strandedness": "ds",
            "stick_ends": None}
        else:
            return {
            "sequence": assembly,
            "topology": "linear",
            "strandedness": "ds",
            "sticky_ends": None
            }
    else:
        return {
            "sequence": assembly,
            "topology": "linear",
            "strandedness": "ds",
            "sticky_ends": None
            }

'''
def assemble_fragments(list_of_subseqs,
                       min_overlap_len = 15,
                       max_overlap_len = 30,
                       product_topo = "linear"):
    #Note: presently requires that all fragments be in the same orientation,
    #in order, and all uppercase (or at least matching cases for the annealing portions)
    #Will fail if one or more provided fragments is not part of the final product
    list_of_subseqs = [seq.upper() for seq in list_of_subseqs]
    print(list_of_subseqs)
    if product_topo not in ["linear"]:
        sys.stderr.write("right now this only works for linear assemblies. Should add circular framework\n")
        sys.exi()
    junctions = []
    for i in range(len(list_of_subseqs)):
        for ii in range(len(list_of_subseqs)):
            for j in range(min_overlap_len, max_overlap_len+1):
                #for seq1, cut slice of j bases from the back end
                #for seq2, cut slice of j bases from the front end
                #if slices match, stich sequences together
                if list_of_subseqs[i][-j:] == list_of_subseqs[ii][:j]:
                    junctions.append((i, ii, j, list_of_subseqs[i][-j:].upper()))
    if junctions == []:
        sys.stderr.write("Unable to find overlap between fragments during assembly, exiting...\n")
        sys.exit()
    junctions = sorted(junctions)
    final_product = ['|']
    for j in junctions:
        #sys.stderr.write('\t'.join([str(k) for k in j]) + '\n')
        if j[1] - j[0] != 1:
            sys.stderr.write("Seq fragments appear to be out of order, please rearrange\n")
            sys.exit()
        if list_of_subseqs[j[0]][:-j[2]].startswith(final_product[-1]):
            final_product.append(list_of_subseqs[j[0]][len(final_product[-1]):-j[2]])
        else:
            final_product.append(list_of_subseqs[j[0]][:-j[2]])
        final_product.append(j[3])
        #sys.stderr.write(''.join(final_product) + '\n')
    final_product.append(list_of_subseqs[junctions[-1][1]][junctions[-1][2]:])
    #sys.stderr.write(''.join(final_product) + '\n')
    #print(final_product)
    return ''.join(final_product).strip('|')
'''
translation_table = {    
    'ACC': "T", 'ATG': "M", 'ACA': "T", 
    'ACG': "T", 'ATC': "I", 'AAC': "N", 
    'ATA': "I", 'AGG': "R", 'CCT': "P", 
    'CTC': "L", 'AGC': "S", 'AAG': "K", 
    'AGA': "R", 'CAT': "H", 'AAT': "N", 
    'ATT': "I", 'CTG': "L", 'CTA': "L", 
    'ACT': "T", 'CAC': "H", 'AAA': "K", 
    'CCG': "P", 'AGT': "S", 'CCA': "P", 
    'CAA': "Q", 'CCC': "P", 'TAT': "Y", 
    'GGT': "G", 'TGT': "C", 'CGA': "R", 
    'CAG': "Q", 'CGC': "R", 'GAT': "D", 
    'CGG': "R", 'CTT': "L", 'TGC': "C", 
    'GGG': "G", 'TAG': "*", 'GGA': "G", 
    'TAA': "*", 'GGC': "G", 'TAC': "Y", 
    'GAG': "E", 'TCG': "S", 'TTT': "F", 
    'GAC': "D", 'CGT': "R", 'GAA': "E", 
    'TCA': "S", 'GCA': "A", 'GTA': "V", 
    'GCC': "A", 'GTC': "V", 'GCG': "A", 
    'GTG': "V", 'TTC': "F", 'GTT': "V", 
    'GCT': "A", 'TTA': "L", 'TGA': "*", 
    'TTG': "L", 'TCC': "S", 'TGG': "W", 
    'TCT': "S"}

def translate(DNA_sequence, frame = 0):
    DNA_sequence = DNA_sequence.upper()
    protein_sequence_list = []
    for i in range(frame,len(DNA_sequence),3):
        protein_sequence_list.append(translation_table[DNA_sequence[i:i+3]])
    return ''.join(protein_sequence_list)

codon_synonyms = {
    'CTT': ['CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
    'ATG': [],
    'AAG': ['AAA'],
    'AAA': ['AAG'],
    'ATC': ['ATA', 'ATT'],
    'AAC': ['AAT'],
    'ATA': ['ATC', 'ATT'],
    'AGG': ['AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
    'CCT': ['CCG', 'CCA', 'CCC'],
    'CTC': ['CTG', 'CTA', 'CTT', 'TTA', 'TTG'],
    'AGC': ['AGT', 'TCG', 'TCA', 'TCC', 'TCT'],
    'ACA': ['ACC', 'ACG', 'ACT'],
    'AGA': ['AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'CAT': ['CAC'],
    'AAT': ['AAC'],
    'ATT': ['ATC', 'ATA'],
    'CTG': ['CTA', 'CTC', 'CTT', 'TTA', 'TTG'],
    'CTA': ['CTG', 'CTC', 'CTT', 'TTA', 'TTG'],
    'ACT': ['ACC', 'ACA', 'ACG'],
    'CAC': ['CAT'],
    'ACG': ['ACC', 'ACA', 'ACT'],
    'CAA': ['CAG'],
    'AGT': ['AGC', 'TCG', 'TCA', 'TCC', 'TCT'],
    'CAG': ['CAA'],
    'CCG': ['CCT', 'CCA', 'CCC'],
    'CCC': ['CCT', 'CCG', 'CCA'],
    'TAT': ['TAC'],
    'GGT': ['GGG', 'GGA', 'GGC'],
    'TGT': ['TGC'],
    'CGA': ['AGG', 'AGA', 'CGC', 'CGG', 'CGT'],
    'CCA': ['CCT', 'CCG', 'CCC'],
    'TCT': ['AGC', 'AGT', 'TCG', 'TCA', 'TCC'],
    'GAT': ['GAC'],
    'CGG': ['AGG', 'AGA', 'CGA', 'CGC', 'CGT'],
    'TTT': ['TTC'],
    'TGC': ['TGT'],
    'GGG': ['GGT', 'GGA', 'GGC'],
    'TAG': ['TAA', 'TGA'],
    'GGA': ['GGT', 'GGG', 'GGC'],
    'TAA': ['TAG', 'TGA'],
    'GGC': ['GGT', 'GGG', 'GGA'],
    'TAC': ['TAT'],
    'GAG': ['GAA'],
    'TCG': ['AGC', 'AGT', 'TCA', 'TCC', 'TCT'],
    'TTA': ['CTG', 'CTA', 'CTC', 'CTT', 'TTG'],
    'GAC': ['GAT'],
    'TCC': ['AGC', 'AGT', 'TCG', 'TCA', 'TCT'],
    'GAA': ['GAG'],
    'TCA': ['AGC', 'AGT', 'TCG', 'TCC', 'TCT'],
    'GCA': ['GCC', 'GCG', 'GCT'],
    'GTA': ['GTC', 'GTG', 'GTT'],
    'GCC': ['GCA', 'GCG', 'GCT'],
    'GTC': ['GTA', 'GTG', 'GTT'],
    'GCG': ['GCA', 'GCC', 'GCT'],
    'GTG': ['GTA', 'GTC', 'GTT'],
    'TTC': ['TTT'],
    'GTT': ['GTA', 'GTC', 'GTG'],
    'GCT': ['GCA', 'GCC', 'GCG'],
    'ACC': ['ACA', 'ACG', 'ACT'],
    'TGA': ['TAG', 'TAA'],
    'TTG': ['CTG', 'CTA', 'CTC', 'CTT', 'TTA'],
    'CGT': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG'],
    'TGG': [],
    'CGC': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT']}