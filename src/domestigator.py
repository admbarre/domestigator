from typing import List, Tuple
from itertools import product
from itertools import chain
import pandas as pd
import re
from Bio.Seq import Seq
import math

DEBUG = True # so it's come to this


class WildSite:
    def __init__(self,seq,site,index):
        self.seq = seq
        self.site = site
        self.index = index
        self.orf_seq = self.get_orf_seq()

    def get_orf_seq(self):
        cut_len = len(self.site)

        orf_index = self.index % 3 # Tells us the offset from the original ORF (0,1,2, 3=0)
        start = max(self.index - orf_index,0) # clamp value at start of the seq
        end = self.index + cut_len
        adjusted_end = math.ceil(end/3)*3 # rounds up to the next codon triplet

        end = min(len(self.seq),    # clamps the end to the end of the seq
                  adjusted_end)
        orf_seq = self.seq[start:end]
        return orf_seq

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
        
    # NOTE: this is the old section when I was approaching
    # with a greedy method. Idk if this is still the way
    def possible_seqs(self,seq: str) -> List[str]:
        residues = Seq(seq).translate()
        codons = [self.res2cod[res] for res in residues]
        seqs = ["".join(p) for p in product(*codons)]
        non_cutters = [seq for seq in seqs if len(self.find_gg(seq)) == 0]
        return non_cutters
    def score_seq(self,seq: str) -> float:
        codons = self.codon_list(seq)
        scores = [1 if codon in self.human_opt.values() else 0 for codon in codons ]
        score = sum(scores)/len(scores)
        return score
    def get_optimal(self,seq_list: List[str]) -> str:
        scored = reversed(sorted([(seq,self.score_seq(seq)) for seq in seq_list],key = lambda x: x[1]))
        top = max(scored,key=lambda x: x[1])[0]
        return top
