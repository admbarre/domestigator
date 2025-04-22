from typing import List, Tuple
from itertools import product
from itertools import chain
import pandas as pd
import re
from Bio.Seq import Seq
import math
import Levenshtein

DEBUG = True # so it's come to this

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
    def find_gg(self,seq):
        cuts = [enz.find_cuts(seq) for enz in self.enzymes]
        cuts = list(chain(*cuts))
        return cuts
    def has_gg(self,seq):
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
        

    def get_orf_seq(self):
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
    def get_checks(self):
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

    def get_seq_freq(self,codons):
        '''
        Converts a list of codons with their frequencies into a seq with a 
        combined frequency.
        '''
        codons,freqs = zip(*codons)
        freq_prod = math.prod(freqs)
        seq = "".join(codons)
        return (seq,freq_prod)


    def possible_seqs(self) -> List[str]:
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
    def domesticate(self,seq):
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

class Optimizer:
    # Lists are ordered in terms of most frequent
    human_optimized_codons = {
        "A": ["GCG", "GCA", "GCC", "GCT"],  # Alanine
        "C": ["TGC", "TGT"],  # Cysteine
        "D": ["GAC", "GAT"],  # Aspartic Acid
        "E": ["GAG", "GAA"],  # Glutamic Acid
        "F": ["TTC", "TTT"],  # Phenylalanine
        "G": ["GGC", "GGA", "GGG", "GGT"],  # Glycine
        "H": ["CAC", "CAT"],  # Histidine
        "I": ["ATC", "ATT", "ATA"],  # Isoleucine
        "K": ["AAG", "AAA"],  # Lysine
        "L": ["CTG", "CTC", "TTG", "CTA", "TTA", "CTT"],  # Leucine
        "M": ["ATG"],  # Methionine
        "N": ["AAC", "AAT"],  # Asparagine
        "P": ["CCG", "CCA", "CCC", "CCT"],  # Proline
        "Q": ["CAG", "CAA"],  # Glutamine
        "R": ["CGC", "CGG", "CGA", "CGT", "AGA", "AGG"],  # Arginine
        "S": ["TCG", "TCC", "AGC", "AGT", "TCA", "TCT"],  # Serine
        "T": ["ACC", "ACA", "ACG", "ACT"],  # Threonine
        "V": ["GTG", "GTC", "GTA", "GTT"],  # Valine
        "W": ["TGG"],  # Tryptophan
        "Y": ["TAC", "TAT"],  # Tyrosine
    }

    human_opt = {
        'A': 'GCC',  # Alanine
        'C': 'TGC',  # Cysteine
        'D': 'GAC',  # Aspartic Acid
        'E': 'GAG',  # Glutamic Acid
        'F': 'TTC',  # Phenylalanine
        'G': 'GGC',  # Glycine
        'H': 'CAC',  # Histidine
        'I': 'ATC',  # Isoleucine
        'K': 'AAG',  # Lysine
        'L': 'CTG',  # Leucine
        'M': 'ATG',  # Methionine (Start)
        'N': 'AAC',  # Asparagine
        'P': 'CCC',  # Proline
        'Q': 'CAG',  # Glutamine
        'R': 'CGC',  # Arginine
        'S': 'AGC',  # Serine
        'T': 'ACC',  # Threonine
        'V': 'GTC',  # Valine
        'W': 'TGG',  # Tryptophan
        'Y': 'TAC',  # Tyrosine
        '*': 'TGA'   # Stop Codon
    }

    res2cod = {
        'F': ['TTT', 'TTC'],  # Phenylalanine
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine
        'M': ['ATG'],  # Methionine (Start)
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine
        'Y': ['TAT', 'TAC'],  # Tyrosine
        '*': ['TAA', 'TAG', 'TGA'],  # Stop codons
        'H': ['CAT', 'CAC'],  # Histidine
        'Q': ['CAA', 'CAG'],  # Glutamine
        'N': ['AAT', 'AAC'],  # Asparagine
        'K': ['AAA', 'AAG'],  # Lysine
        'D': ['GAT', 'GAC'],  # Aspartic Acid
        'E': ['GAA', 'GAG'],  # Glutamic Acid
        'C': ['TGT', 'TGC'],  # Cysteine
        'W': ['TGG'],  # Tryptophan
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
        'G': ['GGT', 'GGC', 'GGA', 'GGG']  # Glycine
    }

    def __init__(self):
        self.domesticated = self.domesticate_cuts()

    def new_find_cuts(self,seq) -> List:
        '''
        Takes a DNA sequence and returns a list containing cut sites.
        '''
    def domesticate_cuts(self,seq):
        '''
        Walks through each cut found for the sequence and attempts to break the
        cut without introducing a new one by optimizing first, and then
        last resort deoptimizing.
        '''
        cuts = self.find_gg(seq)
        domestications = []
        domesticated = None
        for cut,index in cuts:
            if DEBUG:
                print(f"Cut: {cut}")
            checks,subseq,new_index = self.get_orf_subseq(seq,cut,index)
            domesticated = self.opt_domesticate(subseq,checks)
            if not domesticated:
                domesticated = self.deopt_domesticate(subseq,checks)
                if not domesticated:
                    print("could not fix even by deopt...")
                    print("---"*8)
                    domesticated = None
            print("---"*8)
            print()
            dom = (cut,index,subseq,new_index,domesticated)
            domestications.append(dom)
        for dom in domestications:
            cut,index,subseq,new_index,change = dom

            # NOTE: I don't know if this guarantees that we will find this specific
            # subsequence and change it, it is possible this subseq exists elswhere
            # that isn't part of a cute site
            if change == None:
                return None
            domesticated = seq.replace(subseq,change)
        return domesticated

    def find_cuts(self,enz,seq):
        '''
        Finds all cut sites for a given restriction site in the provided
        orientation.
        '''
        cuts = [(enz, m.start()) for m in re.finditer(enz, seq)]
        return cuts
        
    def find_gg(self,seq):
        '''
        Checks the three TypeIIs enzymes we are using in their
        forward and reverse strand orientations. Returns a list
        of these cuts and their indices in their original seq.
        '''
        # We are checking the rev comp of enzyme so we don't have to rev comp the seq and mess up indexing
        bsmbi = "CGTCTC"
        bsmbi_rev = "GCAGAG"
        
        bsai = "GGTCTC"
        bsai_rev = "GAGACC"
        
        paqci = "CACCTGC"
        paqci_rev = "GCAGGTG"

        # TODO: PAQCI IS CURRENTLY BROKEN WITH GETTING ORFs
        enzymes = [bsmbi,bsmbi_rev,bsai,bsai_rev,paqci,paqci_rev]
        #enzymes = [bsmbi,bsmbi_rev,bsai,bsai_rev]
        
        gg_cuts = [self.find_cuts(enz,seq) for enz in enzymes]
        
        return list(chain(*gg_cuts))
        
    def codon_list(self,seq: str) -> List[str]:
        '''
        Returns a string sequence of DNA as a list of string codons.
        '''
        codons = [seq[i:i+3] for i in range(0,len(seq),3)]
        return codons
        
    # TODO: I think it would be better to separate out the checking
    # portion from this method. Ideally we would just get the subseq
    # that we can manipulate as well as the index where it starts
    def get_orf_subseq(self,seq: str,cut: str,index: int):
        '''
        Given a cut site with cut and index and the original sequence,
        provides the cut site in the correct ORF as well as neighboring sequences
        to ensure that an edit doesn't introduce a new cut site by accident.
        '''
        cut_len = math.ceil(len(cut)/3)*3
        orf_index = index % 3 # Tells us the offset from the original ORF (0,1,2, 3=0)
        start = max(index - orf_index,0) # clamp value at 0

        # NOTE: This breaks if our restriction site isn't a multiple of 3!
        # TODO: figure out how to generalize to allow PaqcI or just do those manually
        end = min(start + cut_len ,len(seq)) # clamp value at length of seq
        
        if len(cut) == 7 and DEBUG: # PaqcI
            print("PaqcI")
            print(f"Start: {start} Index: {index} ORF: {orf_index} End: {end}")
        # How far left/right of the subseq we need to test to not accidentally create a new site
        margin = math.ceil(len(cut)/3)*3
        left = max(start - margin,0)
        right = min(end + margin,len(seq))
        
        orf_cut = seq[start:end]
        check_left = seq[left:start]
        check_right = seq[end:right]
        checks = (check_left,check_right)
        if DEBUG:
            print(check_left,orf_cut,check_right)
        return checks,orf_cut,start

    def deopt_domesticate(self,subseq: str,check: Tuple) -> str:
        '''
        This takes a sequence and attemps to break cut site without introducing a new one
        by changing one codon at a time. Only occurs after opt_domesticate has failed and
        we resort to deoptimizing codons.
        '''
        print("vvv"*8)
        print(f"- Attempting to fix through deopt")
        left, right = check
        residues = Seq(subseq).translate()
        original_codons = self.codon_list(subseq)
        
        # Iterating and trying to optimize one codon at a time until site is domesticated
        for i,original in enumerate(original_codons):        
            print(f"Testing codon position: {i}")
            aa = residues[i]

            # List in order of optimization, skip top since it's already opt
            codons_to_test = self.human_optimized_codons[aa][1:]
            for codon in codons_to_test:
                print(f"> {codon}")
                # shallow copy and update
                new_codons = original_codons[:]
                new_codons[i] = codon
                new_seq = "".join(new_codons)

                print(self.codon_list(new_seq))
                # Append left and right neighbors to check
                check_seq = left + new_seq + right
                new_cuts = len(self.find_gg(new_seq))
                if new_cuts == 0:
                    print("Replacement found!")
                    return new_seq
        return None
        
    def opt_domesticate(self,subseq: str, check: Tuple) -> str:
        '''
        This takes a sequence and attemps to break cut site without introducing a new one
        by changing one codon at a time to an optimal codon.
        '''
        print("- Attempting to fix through optimizing")
        left, right = check
        residues = Seq(subseq).translate()
        original_codons = self.codon_list(subseq)
        
        if DEBUG:
            print(f"Old seq: {subseq}")
            print(original_codons)

        # Iterating and trying to optimize one codon at a time until site is domesticated
        for i,original in enumerate(original_codons):        
            # Is possible codon is already opt, move on to next
            if original in self.human_opt.values():
                continue

            # else we optimize the codon and check to see if site is destroyed and no new site created
            aa = residues[i]
            opt_codon = self.human_opt[aa]

            # shallow copy and update
            new_codons = original_codons[:]
            new_codons[i] = opt_codon
            new_seq = "".join(new_codons)
            print(self.codon_list(new_seq))

            # Append left and right neighbors to check
            check_seq = left + new_seq + right
            new_cuts = len(self.find_gg(new_seq))

            # Inspect
            if DEBUG:
                print(f"New seq: {check_seq}")
                print(self.codon_list(new_seq))
                print(f"New cuts generated: {new_cuts}")
            if new_cuts == 0:
                print("Replacement found!")
                return new_seq
            if DEBUG:
                print("Cut found; trying again")
        print("Codons already optimized, need to deopt to fix.")
        print("---"*8)
        return None
        
    def score_seq(self,seq: str) -> float:
        codons = self.codon_list(seq)
        scores = [1 if codon in self.human_opt.values() else 0 for codon in codons ]
        score = sum(scores)/len(scores)
        return score
    def get_optimal(self,seq_list: List[str]) -> str:
        scored = reversed(sorted([(seq,self.score_seq(seq)) for seq in seq_list],key = lambda x: x[1]))
        top = max(scored,key=lambda x: x[1])[0]
        return top
