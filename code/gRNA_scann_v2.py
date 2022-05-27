#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
from collections import Counter
import argparse
import os
from Bio.Blast.Applications import NcbiblastpCommandline
import distance
import scann
from Bio.SeqUtils import GC
import doench_predict

#Functions
def read_fasta(name):
    fasta_seqs = SeqIO.parse(open(name),'fasta')
    data = []
    for fasta in fasta_seqs:
        if not ('mitochondrion' in fasta.description or 'plastid' in fasta.description or 'chloroplast' in fasta.description): #Considering only nuclear chromosomes, removing mitochondrial and plastid/chloroplast genomes
            data.append([fasta.id, str(fasta.seq).strip().upper()])
            
    return data
    
def pam_to_search(pam, iupac_code):
    pam_string = list(pam)
    pam_seq = []
    for i in range(len(pam_string)):
        curr_char = iupac_code[pam_string[i]].split('|')
        if i == 0:
            pam_seq = curr_char
        else:
            curr_pam_seq = []
            for j in range(len(curr_char)):
                for k in range(len(pam_seq)):
                    curr_pam_seq.append(pam_seq[k]+curr_char[j])
                    
            pam_seq = curr_pam_seq
    
    return pam_seq

def one_hot(seq_list):
    d_len = len(seq_list[0])
    one_hot_list = np.zeros((len(seq_list), 4*d_len))
    mapping = dict(zip("ACGT", range(4))) 
    for i in range(len(seq_list)):
        one_hot_list[i] = np.squeeze(np.eye(4)[[mapping[i] for i in seq_list[i]]].reshape(-1,1))

    return np.array(one_hot_list).astype('float32')

def check_intergenic(gtable, chrom, loc, strand, chrom_len, intergenic_space):
    flag = True
    
    #getting the gene list present on the chromsomose of interest
    gloi = gtable.loc[gtable['Accession'] == chrom].reset_index(drop=True)
    
    if strand == '+':
        for i in range(len(gloi)):
            if loc > gloi['Start'][i] - intergenic_space and loc < gloi['Stop'][i] + intergenic_space:
                flag = False
                break
        
    else:
        for i in range(len(gloi)):
            if chrom_len - loc > gloi['Start'][i] - intergenic_space and chrom_len - loc < gloi['Stop'][i] + intergenic_space:
                flag = False
                break
    
    return flag

def get_gene_info(gtable, chrom, loc, strand, chrom_len, gdenslen):
    #getting the gene list present on the chromsomose of interest
    gloi = gtable.loc[gtable['Accession'] == chrom].sort_values('Start').reset_index(drop=True)
    if len(gloi) > 0:
        if strand == '+':
            curr_loc = loc
        else:
            curr_loc = chrom_len - loc

        #left gene, right gene is defined with respect to the top strand [gRNA orientation not taken into account]
        if curr_loc < gloi['Start'][0]:
            right_gene = gloi['Locus tag'][0]
            intg_region_size = '-'
            relative_orient = '-'
            left_gene = '-'

        elif curr_loc > list(gloi['Stop'])[-1]:
            left_gene = list(gloi['Locus tag'])[-1]
            intg_region_size = '-'
            relative_orient = '-'
            right_gene = '-'

        else:
            stop_index = len([x for x in list(gloi['Stop']) if x < curr_loc]) - 1
            intg_region_size = gloi['Start'][stop_index+1] - gloi['Stop'][stop_index]
            left_gene = gloi['Locus tag'][stop_index]
            right_gene = gloi['Locus tag'][stop_index+1]
            
            left_gene_strand = gloi['Strand'][stop_index]
            right_gene_strand = gloi['Strand'][stop_index+1]
            
            if left_gene_strand == right_gene_strand:
                relative_orient = 'tandem'
            elif left_gene_strand == '+' and right_gene_strand == '-':
                relative_orient = 'convergent'
            else:
                relative_orient = 'divergent'

        if curr_loc - gdenslen < 0:
            l_limit = 0
        else: 
            l_limit = curr_loc - gdenslen
            
        if curr_loc + gdenslen > chrom_len:
            u_limit = chrom_len
        else:
            u_limit = curr_loc + gdenslen
            
        gene_density = len(list(set([i for i,v in enumerate(list(gloi['Start'])) if v > l_limit]) & set([i for i,v in enumerate(list(gloi['Stop'])) if v < u_limit])))
        
    else:
        intg_region_size = '-'
        left_gene = '-'
        right_gene = '-'
        relative_orient = '-'
        gene_density = '-'
        
    return [intg_region_size, left_gene, right_gene, relative_orient, gene_density]    

def grna_search(genome, pam_l, glen, orient):
    grna_list = []
    #'for' loop for pam sequences
    for i in range(len(pam_l)):
        curr_pam = pam_l[i]
        
        #'for' loop for chromosomes
        for j in range(len(genome)):
            #top_strand
            chrom_seq = genome[j][1]
            if orient == '3prime':
                curr_grna_loc = [x - glen for x in [m.start() for m in re.finditer(curr_pam, chrom_seq)]]
            else: 
                curr_grna_loc = [m.start() for m in re.finditer(curr_pam, chrom_seq)]
                
            for k in range(len(curr_grna_loc)):
                if curr_grna_loc[k] > -1 and curr_grna_loc[k] < len(chrom_seq) - glen - len(curr_pam):
                    grna_list.append([chrom_seq[curr_grna_loc[k]:curr_grna_loc[k] + glen + len(curr_pam)], genome[j][0], curr_grna_loc[k], '+', len(chrom_seq)])
                    
            #bottom_strand
            chrom_seq = str(Seq(genome[j][1]).reverse_complement())
            if orient == '3prime':
                curr_grna_loc = [x - glen for x in [m.start() for m in re.finditer(curr_pam, chrom_seq)]]
            else: 
                curr_grna_loc = [m.start() for m in re.finditer(curr_pam, chrom_seq)]
                
            for k in range(len(curr_grna_loc)):
                if curr_grna_loc[k] > -1 and curr_grna_loc[k] < len(chrom_seq) - glen - len(curr_pam):
                    grna_list.append([chrom_seq[curr_grna_loc[k]:curr_grna_loc[k] + glen + len(curr_pam)], genome[j][0], curr_grna_loc[k], '-', len(chrom_seq)])
    
    return grna_list

def grna_filter(grna_list, glen, pam, orient, seedlen, re_grna_list, polyG_len, polyT_len, edit_dist, gtable, intergenic_space, gdenslen, ambiguous_nucleotides):
    
    #get grna occurrence table (without PAM)
    if orient == '3prime':
        grna_without_pam = [item[0][0:glen] for item in grna_list]    
    elif orient == '5prime':
        grna_without_pam = [item[0][len(pam):] for item in grna_list]  
        
    #remove ambiguous nucleotides
    grna_without_pam_f = [word for word in grna_without_pam if not any(bad in word for bad in ambiguous_nucleotides)]
    
    grna_occurrence = pd.DataFrame.from_dict(Counter(grna_without_pam_f), orient='index').reset_index()
    grna_occurrence.columns = ['grna', 'frequency']

    #get all guide sequences occuring in genome with duplicates removed (definition of duplicate: sequence occurring in multiple places)
    complete_grna_library_wo_pam = list(grna_occurrence['grna']) 
    
    #get all seed sequences occuring in genome
    if orient == '3prime':
        seed_region = [item[glen-seedlen:] for item in complete_grna_library_wo_pam]    
    elif orient == '5prime':
        seed_region = [item[:seedlen] for item in complete_grna_library_wo_pam]
    
    seed_occurrence = pd.DataFrame.from_dict(Counter(seed_region), orient='index').reset_index()
    seed_occurrence.columns = ['seed_seq', 'frequency']
    complete_seed_library_with_unique_seed = list(seed_occurrence.loc[seed_occurrence['frequency'] == 1, 'seed_seq'])
    
    #get all guide sequences occuring in genome once
    grna_once_wo_pam_library_f = list(grna_occurrence.loc[grna_occurrence['frequency'] == 1]['grna'])

    #gRNA RE check
    re_to_check_grna_pre = re_grna_list.split(',')   
    
    if not '' in re_to_check_grna_pre:
        #Taking care of ambiguous nucleotides
        re_to_check_grna = []
        for i in range(len(re_to_check_grna_pre)):
            re_to_check_grna.extend(pam_to_search(re_to_check_grna_pre[i] ,iupac_code))
        
        #checking provided recognition sequence
        grna_wo_pam_f = [w for w in grna_once_wo_pam_library_f if not any(b in w for b in re_to_check_grna)] #f stands for filter
        
        #checking reverse complement of provided recognition sequence
        revcom_re_to_check_grna = []
        for i in range(len(re_to_check_grna)):
            revcom_re_to_check_grna.append(str(Seq(re_to_check_grna[i]).reverse_complement()))
            
        grna_wo_pam_f = [w for w in grna_wo_pam_f if not any(b in w for b in revcom_re_to_check_grna)]
    else: 
        grna_wo_pam_f = grna_once_wo_pam_library_f
        
    #gRNA polyG check
    if not polyG_len == 0:
        polyG_to_check = "G" * polyG_len
        grna_wo_pam_f = [w for w in grna_wo_pam_f if polyG_to_check not in w]
        
        polyC_to_check = "C" * polyG_len
        grna_wo_pam_f = [w for w in grna_wo_pam_f if polyC_to_check not in w]
        
    #gRNA polyT check
    if not polyT_len == 0:
        polyT_to_check = "T" * polyT_len
        grna_wo_pam_f = [w for w in grna_wo_pam_f if polyT_to_check not in w]
        
        polyA_to_check = "A" * polyT_len
        grna_wo_pam_f = [w for w in grna_wo_pam_f if polyA_to_check not in w]
        
    #get guides sequences with a unique seed
    if orient == '3prime':
        grna_wo_pam_f_seed_region = [item[glen-seedlen:] for item in grna_wo_pam_f]    
    elif orient == '5prime':
        grna_wo_pam_f_seed_region = [item[:seedlen] for item in grna_wo_pam_f]    

    unique_seed_library = list(set(grna_wo_pam_f_seed_region) & set(complete_seed_library_with_unique_seed))
    seed_to_compare_dict = dict(zip(grna_wo_pam_f_seed_region,range(0,len(grna_wo_pam_f_seed_region))))
    grna_wo_pam_f_index = [seed_to_compare_dict.get(key) for key in unique_seed_library]
    grna_wo_pam_us_f = [grna_wo_pam_f[i] for i in grna_wo_pam_f_index]
    
    #off-target using unique_seed_library and complete_seed_library
    xq = one_hot(grna_wo_pam_us_f)
    xb = one_hot(complete_grna_library_wo_pam)
    norm_xb = xb/np.linalg.norm(xb, axis=1)[:, np.newaxis]
    
    searcher = scann.scann_ops_pybind.builder(norm_xb, 10, "dot_product").tree(
    num_leaves=2000, num_leaves_to_search=100, training_sample_size=5000000).score_ah(
    2, anisotropic_quantization_threshold=0.2).reorder(100).build()
    
    neighbors, distances = searcher.search_batched(xq, leaves_to_search=250, pre_reorder_num_neighbors=250)
    
    unique_grna_library_mm = []
    k_to_check = 3 #knearest neighbor search
    for i in range(len(xq)):
        knn_dist = []
        for j in range(k_to_check):
            knn_dist.append(distance.hamming(grna_wo_pam_us_f[i], complete_grna_library_wo_pam[neighbors[i][j+1]]))
    
        if all(i > edit_dist - 1 for i in knn_dist):
            unique_grna_library_mm.append(grna_wo_pam_us_f[i])
    
    if orient == '3prime':
        grna_to_compare = [item[0][:glen] for item in grna_list]    
    elif orient == '5prime':
        grna_to_compare = [item[0][len(pam):len(pam)+glen] for item in grna_list]    

    grna_to_compare_dict = dict(zip(grna_to_compare,range(0,len(grna_to_compare))))
    grna_index = [grna_to_compare_dict.get(key) for key in unique_grna_library_mm]
    grna_list_mm = [grna_list[i] for i in grna_index]
    
    #select intergenic gRNA
    if len(grna_list_mm) > 0:
        index_to_keep = []
        for i in range(len(grna_list_mm)):
            if check_intergenic(gtable, grna_list_mm[i][1], grna_list_mm[i][2], grna_list_mm[i][3], grna_list_mm[i][4], intergenic_space):
                index_to_keep.append(i)

        grna_list_mm_intg = [grna_list_mm[i] for i in index_to_keep]

        #intergenic region size, adjacent genes and gene density
        for i in range(len(grna_list_mm_intg)):
            grna_list_mm_intg[i].extend(get_gene_info(gtable, grna_list_mm_intg[i][1], grna_list_mm_intg[i][2], grna_list_mm_intg[i][3], grna_list_mm_intg[i][4], gdenslen))
    
    return grna_list_mm_intg

def hr_filter(data, glen, pam, genome, hr_len, RE_hr, polyG_hr, polyT_hr):
    for i in range(len(data)):
        
        #get the chromosome seq based on data[i][1]; use genome
        chrom_seq = genome[[chrom_name[0] for chrom_name in genome].index(data[i][1])][1]
        
        #pam length 
        pam_len = len(pam)
        if data[i][3] == '+':
            curr_chrom_seq = chrom_seq
        else:
            curr_chrom_seq = str(Seq(chrom_seq).reverse_complement())
            
        if data[i][2] - hr_len > 0:
            LHR = curr_chrom_seq[data[i][2] - hr_len:data[i][2]]
        else:
            LHR = '-'
            
        if data[i][2] + glen + pam_len + hr_len < data[i][4]:
            RHR = curr_chrom_seq[data[i][2] + glen + pam_len:data[i][2] + glen + pam_len + hr_len]
        else:
            RHR = '-' 
            
        data[i].extend([LHR, RHR])
    
    LHR_list = [lhr_seq[-2] for lhr_seq in data]
    RHR_list = [rhr_seq[-1] for rhr_seq in data]
    
    index_to_remove = []
    #RE_check
    re_to_check_hr_pre = RE_hr.split(',')   
    if not '' in re_to_check_hr_pre:
        #Taking care of ambiguous nucleotides
        re_to_check_hr = []
        for i in range(len(re_to_check_hr_pre)):
            re_to_check_hr.extend(pam_to_search(re_to_check_hr_pre[i] ,iupac_code))
        
        for i in range(len(data)):
            for j in range(len(re_to_check_hr)):
                if re_to_check_hr[j] in LHR_list[i] or re_to_check_hr[j] in RHR_list[i]:
                    index_to_remove.append(i)

    #PolyG_check
    if not polyG_hr == 0:
        polyG_to_check = "G" * polyG_hr
        polyC_to_check = "C" * polyG_hr
        for i in range(len(data)):
            if polyG_to_check in LHR_list[i] or polyC_to_check in RHR_list[i]:
                index_to_remove.append(i)
        
    #polyT check
    if not polyT_len == 0:
        polyT_to_check = "T" * polyT_hr
        polyA_to_check = "A" * polyT_hr
        for i in range(len(data)):
            if polyT_to_check in LHR_list[i] or polyA_to_check in RHR_list[i]:
                index_to_remove.append(i)
    
    gh_data = [i for j, i in enumerate(data) if j not in np.unique(index_to_remove)]
    
    return gh_data

def write_fasta(name, sequence_df):
    out_file = open(name, "w")
    for i in range(len(sequence_df)):
        out_file.write('>' + sequence_df['DEG_id'][i] + '\n')
        out_file.write(sequence_df['Sequence'][i] + '\n')
    out_file.close()

#IUPAC Code
iupac_code = {
  "A": "A",
  "T": "T",
  "G": "G",
  "C": "C",
  "M": "A|C",
  "R": "A|G",
  "W": "A|T",
  "S": "C|G",
  "Y": "C|T",
  "K": "G|T",
  "V": "A|C|G",
  "H": "A|C|T",
  "D": "A|G|T",
  "B": "C|G|T",
  "N": "A|C|G|T",
}

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--Genome', help="Genome filename", required=True)
parser.add_argument('-t', '--Gene_table', help="Gene Table filename", required=True)
parser.add_argument('-out', '--Output_file', help="Name of the output file")
parser.add_argument('-p', '--PAM', type=str, default='NGG', help="A short PAM motif to search for, it may use IUPAC ambiguous alphabet. Default: NGG.", required=True)
parser.add_argument('-o', '--Orientation', choices=['3prime', '5prime'], default='3prime', help="PAM position relative to target: 5prime: [PAM][target], 3prime: [target][PAM]. For example, Cas9 is 3prime. Default: '3prime'.")
parser.add_argument('-l', '--Guide_Length', type=int, choices=range(10, 28, 1), metavar="[10-27]", default=20, help="Length of the guide sequence. Default: 20.")
parser.add_argument('-sl','--Seed_Length', type=int, choices=range(0, 28, 1), metavar="[0-27]", default=10, help='Length of a seed region near the PAM site required to be unique. Default: 10.')
parser.add_argument('--RE_grna', type=str, default='', help='Undesired recognition sequence of restriction enzymes in guide RNA')
parser.add_argument('--polyG_grna', type=int, choices=range(0, 11, 1), metavar="[0-10]", default=0, help='Length of G/C repeats not allowed in the guide sequence. Default value of 0 implies poly_G rule is not applied.')
parser.add_argument('--polyT_grna', type=int, choices=range(0, 11, 1), metavar="[0-10]", default=0, help='Length of T/A repeats not allowed in the guide sequence. Default value of 0 implies poly_T rule is not applied.')
parser.add_argument('--intspace', type=int, default=300, help='Minimum distance of gRNA from any gene. Default is 300bp.')
parser.add_argument('--edit_dist', type=int, default=6, choices=range(0, 11, 1), metavar="[0-10]",  help='Minimum number of mismatch in the guide region to classify a guide sequence as candidate. Default value is 6.')
parser.add_argument('--gene_density_len', type=int, default=10000, help='Size of the region from the gRNA site to calculate gene density. Default is 10000bp.')
parser.add_argument('-hr_l', '--HR_Length', type=int, choices=range(5, 201, 1), metavar="[5-500]", default=50, help="Length of the Homology Arm. Default: 50bp.")
parser.add_argument('--RE_hr', type=str, default='', help='Undesired recognition sequence of restriction enzymes in the homology arm.')
parser.add_argument('--polyG_hr', type=int, choices=range(0, 11, 1), metavar="[0-10]", default=0, help='Length of G/C repeats not allowed in the homology arm. Default value of 0 implies poly_G rule is not applied.')
parser.add_argument('--polyT_hr', type=int, choices=range(0, 11, 1), metavar="[0-10]", default=0, help='Length of T/A repeats not allowed in the homology arm. Default value of 0 implies poly_T rule is not applied.')
parser.add_argument('--backbone_complementarity_check', type=str, default='', help='Check if guide RNA will form secondary structure with the backbone.')
parser.add_argument('--protein_file', type=str, default='', help="Fasta file containing protein sequences")
parser.add_argument('--blast_org', type=str, default='',  help="Organisms to blast proteins against to identify probable essential genes.")

args = parser.parse_args()
genome_file = args.Genome
gene_file = args.Gene_table
output_file = args.Output_file
pam = args.PAM
orient = args.Orientation
glen = args.Guide_Length
seedlen = args.Seed_Length
re_grna_list = args.RE_grna
polyG_len = args.polyG_grna
polyT_len = args.polyT_grna
intergenic_space = args.intspace
edit_dist = args.edit_dist
gdenslen = args.gene_density_len
hr_len = args.HR_Length
re_hr_list = args.RE_hr
polyG_hr = args.polyG_hr
polyT_hr = args.polyT_hr
protein_file = args.protein_file
org_ge = args.blast_org
backbone_region = args.backbone_complementarity_check

#Data Processing
path = '/home/aboob2/universal_landing_pad/data/'
genome = read_fasta(path + genome_file)
gene_table = pd.read_csv(path + gene_file)
refined_gene_table = gene_table[['Accession', 'Start', 'Stop', 'Strand', 'Locus tag']]
pam_library = pam_to_search(pam,iupac_code)
ambiguous_nucleotides = list(iupac_code.keys())[4:]

#Obtaining harbors
grna_list = grna_search(genome, pam_library, glen, orient)
grna_data = grna_filter(grna_list, glen, pam, orient, seedlen, re_grna_list, polyG_len, polyT_len, edit_dist, refined_gene_table, intergenic_space, gdenslen, ambiguous_nucleotides)

if len(grna_data) > 0:
    grna_hr_data = hr_filter(grna_data, glen, pam, genome, hr_len, re_hr_list, polyG_hr, polyT_hr)

    #Cleaning and Labeling dataframe
    grna_hr_df = pd.DataFrame(grna_hr_data, columns = ['Guide with PAM', 'Accession', 'Location', 'Strand', 'Chromosome Length', 'Intergenic Size', 'Left Gene', 'Right Gene', 'Relative Orientation', 'Gene Density', 'Left HR', 'Right HR'])

    if orient == '3prime':
        guide_seq = grna_hr_df['Guide with PAM'].str[:glen]
        pam_seq = grna_hr_df['Guide with PAM'].str[glen:]
    elif orient == '5prime':
        guide_seq = grna_hr_df['Guide with PAM'].str[len(pam):]
        pam_seq = grna_hr_df['Guide with PAM'].str[:len(pam)]

    grna_hr_df.insert(loc=0, column='PAM', value=pam_seq)
    grna_hr_df.insert(loc=0, column='Guide Sequence', value=guide_seq)
    del grna_hr_df['Guide with PAM']

    chrom_name_df = gene_table.drop_duplicates('Accession').reset_index(drop=True)[['#Name','Accession']]
    chrom_name_list = []
    for i in range(len(grna_hr_df)):
        chrom_name_list.append(chrom_name_df.loc[chrom_name_df['Accession'] == grna_hr_df['Accession'][i], '#Name'].iloc[0])

    grna_hr_df.insert(loc=3, column='Chromosome', value=chrom_name_list)
    del grna_hr_df['Accession']
    
    guide_gc = []
    for i in range(len(grna_hr_df)):
        guide_gc.append(GC(grna_hr_df['Guide Sequence'][i]))
        
    grna_hr_df.insert(loc = 2, column='GC% Guide', value = guide_gc)
    
    self_comp = []
    stem_len = 4
    for i in range(len(grna_hr_df)):
        fwd = grna_hr_df['Guide Sequence'][i]
        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd)-stem_len-1

        folding = 0
        for j in range(0,len(fwd)-stem_len):
            if GC(fwd[j:j+stem_len]) >= 0.5:
                if fwd[j:j+stem_len] in rvs[0:(L-j)] or any([fwd[j:j+stem_len] in item for item in backbone_region]):
                    folding += 1

        self_comp.append(folding)

    grna_hr_df.insert(loc = 3, column='Self-Complementarity', value = self_comp)
    
    if pam == 'NGG' and orient == '3prime':
        on_target_seq = []
        for i in range(len(grna_hr_df)):
            if len(grna_hr_df['Guide Sequence'][i]) < 23:
                on_target_seq.append(grna_hr_df['Right HR'][i][len(grna_hr_df['Guide Sequence'][i])-24:] + grna_hr_df['Guide Sequence'][i] + grna_hr_df['PAM'][i] + grna_hr_df['Right HR'][i][0:3])
            else:
                on_target_seq.append(grna_hr_df['Guide Sequence'][i][0:24] + grna_hr_df['PAM'][i] + grna_hr_df['Right HR'][i][0:3])

        grna_hr_df['On-target Score'] = doench_predict.predict(np.array(on_target_seq), num_threads=1)
    else:
        grna_hr_df['On-target Score'] = 'Not Available'

    if org_ge and protein_file:
        #Adding essentiality information
        proteins_query = path + protein_file
        organism_list = org_ge.split(',')

        #procuring Essential Gene Database file
        deg_file = '/home/aboob2/universal_landing_pad/essential genes/deg.csv'
        deg_database = pd.read_csv(deg_file)

        db = os.path.join(path, protein_file.split('/')[0], 'RefOrg.faa') # BLAST database
        blastout = os.path.join(path, protein_file.split('/')[0],'blast.tab')  # BLAST output

        eg_df = deg_database[deg_database.Organism.isin(organism_list)].reset_index(drop=True).iloc[:,0:2]

        #create RefOrg file 
        ref_org = path + protein_file.split('/')[0] + '/RefOrg.fasta'
        write_fasta(ref_org, eg_df)

        blast_path = '/home/aboob2/miniconda3/pkgs/blast-2.5.0-hc0b0e79_3/bin/'
        #Creating Blast Database
        blastdb_cmd = '{}makeblastdb -in {} -parse_seqids -dbtype prot -out {}'.format(blast_path, ref_org, db)
        os.system(blastdb_cmd)

        #Blast
        cmd_blastp = NcbiblastpCommandline(cmd = blast_path + 'blastp', query = proteins_query, out = blastout, outfmt = 6, db = path + protein_file.split('/')[0] + '/RefOrg.faa')
        stdout, stderr = cmd_blastp()

        results = pd.read_csv(blastout, sep="\t", header=None)
        headers = ['query', 'subject',
                   'pc_identity', 'aln_length', 'mismatches', 'gaps_opened',
                   'query_start', 'query_end', 'subject_start', 'subject_end',
                   'e_value', 'bitscore']

        results.columns = headers
        results_filtered = results.loc[(results['e_value'] < 1e-5) & (results['pc_identity'] >= 50)].reset_index(drop=True)
        eg_loc_df = gene_table[gene_table['Protein product'].isin(np.unique(list(results_filtered['query'])))].reset_index(drop=True)

        chrom_len_array = []
        for i in range(len(chrom_name_df)):
            for j in range(len(genome)):
                if chrom_name_df['Accession'][i] == genome[j][0]:
                    chrom_len_array.append(len(genome[j][1]))

        chrom_name_df['Length'] = chrom_len_array

        chr_eg_zone = []
        for i in range(len(chrom_name_df)):
            curr_chr_eg_data = eg_loc_df.loc[eg_loc_df['#Name'] == chrom_name_df['#Name'][i]].sort_values('Start').reset_index(drop=True)

            zone_info = []
            ini_flag = 0
            if np.shape(curr_chr_eg_data)[0] > 0:
                for j in range(len(curr_chr_eg_data)):
                    if ini_flag == 0:
                        zone_info = str(0) + '-' + str(curr_chr_eg_data['Start'][j])
                        chr_eg_zone.append([chrom_name_df['#Name'][i], zone_info])
                        ini_flag = 1

                    if j == len(curr_chr_eg_data) - 1:
                        zone_info = str(curr_chr_eg_data['Stop'][j]) + '-' + str(chrom_name_df['Length'][i])
                    else:
                        zone_info = str(curr_chr_eg_data['Stop'][j]) + '-' + str(curr_chr_eg_data['Start'][j+1])

                    chr_eg_zone.append([chrom_name_df['#Name'][i], zone_info])
            else:
                zone_info = str(0) + '-' + str(chrom_name_df['Length'][i]) #no essential genes on that chromosome
                chr_eg_zone.append([chrom_name_df['#Name'][i], zone_info])

        chr_eg_zone = pd.DataFrame(chr_eg_zone, columns = ['Chr','Loc'])

        zone = []
        for i in range(len(grna_hr_df)):
            if grna_hr_df['Strand'][i] == '+':
                grna_loc = grna_hr_df['Location'][i]
            else:
                grna_loc = grna_hr_df['Chromosome Length'][i] - grna_hr_df['Location'][i]

            for j in range(np.shape(chr_eg_zone)[0]):
                if chr_eg_zone['Chr'][j] == grna_hr_df['Chromosome'][i]:
                    curr_zone_bound = chr_eg_zone['Loc'][j].split('-')
                    if grna_loc > int(curr_zone_bound[0]) and grna_loc < int(curr_zone_bound[1]):
                           zone.append(j+1)

        grna_hr_es_df = grna_hr_df
        grna_hr_es_df['Zone'] = zone
        pd.DataFrame(grna_hr_es_df).to_csv(path + output_file, index = False) #Safe Harbor Data Output

    else:
        pd.DataFrame(grna_hr_df).to_csv(path + output_file, index = False) #Harbor Data Output
        
else:
    print('No Safe Harbors can be obtained after applying the specified constraints. Try relaxing the edit distance criteria.')
