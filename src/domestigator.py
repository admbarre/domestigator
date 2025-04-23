from itertools import product
from itertools import chain
import pandas as pd
import re
from Bio.Seq import Seq
import math
import Levenshtein

human_codon_usage = {
    'F': [('TTT', 0.46), ('TTC', 0.54)],
    'L': [('TTA', 0.07), ('TTG', 0.13), ('CTT', 0.13), ('CTC', 0.20), ('CTA', 0.07), ('CTG', 0.40)],
    'I': [('ATT', 0.36), ('ATC', 0.47), ('ATA', 0.17)],
    'M': [('ATG', 1.00)],
    'V': [('GTT', 0.18), ('GTC', 0.24), ('GTA', 0.12), ('GTG', 0.46)],
    'S': [('TCT', 0.18), ('TCC', 0.22), ('TCA', 0.15), ('TCG', 0.05), ('AGT', 0.15), ('AGC', 0.25)],
    'P': [('CCT', 0.28), ('CCC', 0.33), ('CCA', 0.27), ('CCG', 0.12)],
    'T': [('ACT', 0.24), ('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.11)],
    'A': [('GCT', 0.27), ('GCC', 0.40), ('GCA', 0.23), ('GCG', 0.11)],
    'Y': [('TAT', 0.43), ('TAC', 0.57)],
    'H': [('CAT', 0.42), ('CAC', 0.58)],
    'Q': [('CAA', 0.27), ('CAG', 0.73)],
    'N': [('AAT', 0.47), ('AAC', 0.53)],
    'K': [('AAA', 0.43), ('AAG', 0.57)],
    'D': [('GAT', 0.46), ('GAC', 0.54)],
    'E': [('GAA', 0.42), ('GAG', 0.58)],
    'C': [('TGT', 0.46), ('TGC', 0.54)],
    'W': [('TGG', 1.00)],
    'R': [('CGT', 0.08), ('CGC', 0.19), ('CGA', 0.11), ('CGG', 0.20), ('AGA', 0.21), ('AGG', 0.21)],
    'G': [('GGT', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)],
    '*': [('TAA', 0.30), ('TAG', 0.24), ('TGA', 0.46)]
}

class REnzyme:
    def __init__(self,name,site):
        self.name = name
        self.site = site
        self.rc_site = str(Seq(site).reverse_complement())
    def find_cuts(self,seq):
        fwd_cuts = [(self.site, m.start()) for m in re.finditer(self.site, seq)]
        rev_cuts = [(self.rc_site, m.start()) for m in re.finditer(self.rc_site, seq)]
        cuts = fwd_cuts + rev_cuts
        return cuts


class GG:
    bsai = REnzyme("BsaI","GGTCTC")
    bsmbi = REnzyme("BsmbI","CGTCTC")
    paqci = REnzyme("PaqcI","CACCTGC")
    enzymes = [bsai,bsmbi,paqci]
    def __init__(self):
        pass
    def find_gg(self,seq: str) -> list:
        cuts = [enz.find_cuts(seq) for enz in self.enzymes]
        cuts = list(chain(*cuts))
        return cuts
    def has_gg(self,seq: str) -> bool:
        cuts = self.find_gg(seq)
        return len(cuts) > 0
        

class WildSite:
    # TODO: rethink this
    # This feels weird
    gg = GG()
    def __init__(self,seq,site,index):
        self.seq = seq
        self.site = site
        self.cut_len = len(site)
        self.index = index
        self.orf_index = index % 3
        self.orf_start = self.index - self.orf_index
        
        self.orf_seq = self.get_orf_seq()
        self.orf_len = len(self.orf_seq)

        self.left = self.right = ""
        if self.orf_index != 0: # if we are in frame we don't need to check flanks
            self.left,self.right = self.get_checks()
        self.residues = Seq(self.orf_seq).translate()
        self.possible = self.possible_seqs()

        # Get the best sequence
        self.optimal = self.possible[0]
        self.domestication = (self.orf_seq,self.optimal,self.orf_start)
        

    def get_orf_seq(self) -> str:
        '''
        Gets ORF containing cut site.
        '''
        start = max(self.index - self.orf_index,0) # clamp value at start of the seq
        end = self.index + self.cut_len
        adjusted_end = math.ceil(end/3)*3 # rounds up to the next codon triplet

        end = min(len(self.seq),    # clamps the end to the end of the seq
                  adjusted_end)
        orf_seq = self.seq[start:end]
        return orf_seq

    # TODO: just realized that the check window should be the length
    # of the longest cutter (PaqcI) in the event that it creates it
    # so we don't miss it
    def get_checks(self) -> tuple[str,str]:
        '''
        Gets the left/right regions to check if edits introduce a new
        cut site. Only len(cut)-1.
        '''
        left_end = self.index - self.orf_index
        left_start = max(left_end - (self.cut_len - 1),0)
        left = self.seq[left_start:left_end]

        right_start = left_end + self.orf_len
        right_end = min(right_start + self.cut_len - 1, len(self.seq))
        right = self.seq[right_start:right_end]
        return left,right

    def get_seq_freq(self,codons: list[str]) -> tuple[str,float]:
        '''
        Converts a list of codons with their frequencies into a seq with a 
        combined frequency.
        '''
        codons,freqs = zip(*codons)
        freq_prod = math.prod(freqs)
        seq = "".join(codons)
        return (seq,freq_prod)


    def possible_seqs(self) -> list[str]:
        '''
        Sorts by possible synonymous sequences by similarity first
        and then by degree of codon optimization.
        '''
        codons = [human_codon_usage[res] for res in self.residues]
        seq_freq = [self.get_seq_freq(p) for p in product(*codons)]
        seq_freq_dist = [(seq,freq,Levenshtein.ratio(seq,self.orf_seq)) for seq,freq in seq_freq]
        ranked = sorted(seq_freq_dist,
                        key=lambda x: (-x[2],x[1]))

        # NOTE: is it safe to drop the weights now?
        ranked_check = [(f"{self.left}{seq[0]}{self.right}",seq[0]) for seq in ranked]
        non_cutters = [seq for check,seq in ranked_check if not self.gg.has_gg(check)]
        diff = len(ranked_check) - len(non_cutters)
        #print(f"Trimmed {diff} cutting sequences")
        return non_cutters

# NOTE: stop coding at 3:34AM my brain is fried
class Gator:
    gg = GG()
    def __init__(self):
        pass
    def domesticate(self,seq: str) -> str:
        cuts = self.gg.find_gg(seq)
        #print(f"Found {len(cuts)} cuts")
        if len(cuts) == 0:
            #print("No cuts found")
            return seq
        wild_sites = [WildSite(seq,cut,index) for cut,index in cuts]
        domestications = [ws.domestication for ws in wild_sites]
        new_seq = seq
        for dom in domestications:
            old,new,index = dom
            #assert len(old) == len(new)
            #print(old,new,index)
            new_seq = new_seq[:index] + new + new_seq[index+len(new):]
        #print("---")
        return new_seq
