# Copyright 2013 Vinothan N. Manoharan and William B. Rogers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Defines NNModel and subclasses, used to calculate thermodynamic parameters
and melting points of nucleic-acid duplexes using the Nearest-Neighbor
model.

Requires Numpy and the Biopython package

.. moduleauthor:: Vinothan N. Manoharan <vnm@seas.harvard.edu>
.. moduleauthor:: William B. Rogers <wrogers@seas.harvard.edu>
"""

import numpy as np
import string
import re
import Bio      # requires Biopython library for sequence processing
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet
import difflib  # for substring matching

def _count_overlapping(subseq, seq):
    """
    Parse a sequence to get the number of overlapping occurrences of a
    subsequence (see http://stackoverflow.com/questions/6844005/). Uses
    python regular expression library.
    """
    return len(re.findall(r'(?=(%s))' % re.escape(subseq), seq))

def _disambiguate(seq):
    """
    Crude utility function to identify a Seq object or string as DNA
    or RNA.  Returns a Seq object with the identified alphabet and all
    letters uppercased.
    """
    useq = seq.upper()
    if isinstance(useq, Seq):
        if ((useq.alphabet == IUPAC.unambiguous_rna) or
            (useq.alphabet == IUPAC.unambiguous_dna)):
            return useq
    if useq.find('U') == -1:
        # assume it's DNA if there's no U
        s = Seq(str(useq), alphabet=IUPAC.unambiguous_dna)
    else:
        s = Seq(str(useq), alphabet=IUPAC.unambiguous_rna)
    return s

def distance_map(seq1, seq2):
    """
    Find differences between seq1 and the complementary part of seq2
    as a function of position along seq2.  Returns the Levenshtein
    distance at each position.  Useful for testing how
    specifically an oligo will attach to a target sequence
    """
    seqmatcher = difflib.SequenceMatcher(None, str(seq1), str(seq2),
                                         autojunk=False)
    

class Duplex:
    """
    Class defining a duplex as a pair of two Biopython Seq objects.
    Duplexes can be DNA-DNA, DNA-RNA, or RNA-RNA.  Mismatches and
    dangling ends are allowed.
    """
    def __init__(self, strand1, strand2=None, align=(0,0)):
        """ 
        Constructs a duplex from two strands.  The duplex contains
        only the subsequences of the two strands that are specified by
        the 'align' parameter. 

        Parameters
        ----------
        strand1 : Biopython Seq object or string
            First strand
        strand2 : Biopython Seq object or string, optional
            Second strand.  If not specified, the reverse
            complement is used.  
        align : tuple of int, optional
            Specifies the beginning and end position of the duplex.
            First number is the position on strand 1 where the duplex
            begins, and the second number is the position on strand 2
            where the duplex ends. Default is (0,0), meaning that the
            5' end of the first strand is the beginning of the duplex,
            and the 5' end of the second strand is the end.

        Notes
        -----
        If a hybrid DNA-RNA duplex is specified, the seq attribute
        will be the DNA strand and the cseq attribute will be the RNA
        strand, regardless of the order in which the strands were
        given.
        """
        s1 = _disambiguate(strand1)
        if _verify_alphabet(s1) is False:
            raise ValueError("Couldn't identify strand 1 %s as "
                             "unambiguous DNA or RNA" %str(strand1))
        if strand2 is None:
            s2 = s1.reverse_complement()
        else:
            s2 = _disambiguate(strand2)
            if _verify_alphabet(s2) is False:
                raise ValueError("Couldn't identify strand 2 %s as "
                                 "as unambiguous DNA or RNA" %str(strand2))

        # calculate dangling ends (TODO: implement)
        i = align[0]
        m = align[1]
        #j = min(len(s1), )
        #n = min(

        # For DNA-RNA hybrids, use the DNA strand as seq and the
        # RNA strand as cseq
        if ((s1.alphabet == IUPAC.unambiguous_rna) and
            (s2.alphabet == IUPAC.unambiguous_dna)):
            self.seq = s2
            self.cseq = s1
        else:
            self.seq = s1
            self.cseq = s2

    def __str__(self):
        """
        Format the two sequences in a duplex in human-readable
        (stacked) form
        """
        duplex_string = "5' {} 3'\n3' {} 5'".format(self.seq,
                                                    self.cseq[::-1]) 
        return duplex_string
        
    def is_self_complementary(self):
        """
        Method to check whether a duplex contains self-complementary
        strands.  Always returns False for DNA-RNA duplexes
        """
        if str(self.seq) == str(self.cseq):
            return True
        else: 
            return False

        
class NNModel:
    """
    Defines common interface and shared code for all variants of the
    nearest-neighbor model (DNA-DNA, RNA-RNA, and DNA-RNA)

    Notes
    -----
    """
    def __init__(self):
        self.Tref = 37+273.15   # reference temperature for thermo data
        self.R = 1.987e-3       # gas constant in kcal/mol/K

        # pairs are written 5' to 3'. pairs is for a DNA strand; rpairs
        # for an RNA strand
        self.pairs = ['AA', 'AT', 'AC', 'AG', 'TA', 'TT', 'TC', 'TG',
                       'CA', 'CT', 'CC', 'CG', 'GA', 'GT', 'GC', 'GG']
        self.rpairs = ['AA', 'AU', 'AC', 'AG', 'UA', 'UU', 'UC', 'UG',
                       'CA', 'CU', 'CC', 'CG', 'GA', 'GU', 'GC', 'GG']

    def calc_thermo(self, duplex, T=37+273.15, saltc=None):
        """
        Virtual function; subclass implementation does
        thermodynamic calculations based on nearest-neighbor model.
        """
        pass
    
    def melting_temp(self, duplex, c=1e-4, saltc=None):
        """
        Returns melting point of duplex in K.  

        Parameters
        ----------
        duplex : Duplex object
            Duplex of interest
        c : float, optional
            total concentration of strand and complement (in mol/L).
        saltc : float, optional
            concentration of sodium ions in mol/L

        Notes
        -----
        Currently assumes strand and complement are in equimolar
        amounts.  The correction factor F is therefore equal to 4 for
        non-self-complementary duplexes.
        """
        if duplex.is_self_complementary():
            sym = 1
        else:
            sym = 4
        dH, dS, dG = self.calc_thermo(duplex, saltc=saltc)
        return dH/(self.R*np.log(c/sym) + dS)
        
        
class DNAModel(NNModel):
    def __init__(self):
        NNModel.__init__(self)

        # Use dictionaries to store the NN data so that there's no
        # confusion about which value belongs to which NN pair.
        # Units of energy are kcal/mole.  Units of entropy are are "e.u",
        # or cal/mol/K.  Need to divide by 1000 before multiplying by T

        # DNA-DNA duplex data.
        # From SantaLucia and Hicks, Annu. Rev. Biophys. Biomol. Struct 2004
        self.dh = {'AA':-7.6, 'AT':-7.2, 'AC':-8.4, 'AG':-7.8,
                   'TA':-7.2, 'TT':-7.6, 'TC':-8.2, 'TG':-8.5,
                   'CA':-8.5, 'CT':-7.8, 'CC':-8.0, 'CG':-10.6,
                   'GA':-8.2, 'GT':-8.4, 'GC':-9.8, 'GG':-8.0}

        self.ds = {'AA':-21.3, 'AT':-20.4, 'AC':-22.4, 'AG':-21.0,
                   'TA':-21.3, 'TT':-21.3, 'TC':-22.2, 'TG':-22.7,
                   'CA':-22.7, 'CT':-21.0, 'CC':-19.9, 'CG':-27.2,
                   'GA':-22.2, 'GT':-22.4, 'GC':-24.4, 'GG':-19.9}

        self.dg = {'AA':-1.00, 'AT':-0.88, 'AC':-1.44, 'AG':-1.28,
                   'TA':-0.58, 'TT':-1.00, 'TC':-1.30, 'TG':-1.45,
                   'CA':-1.45, 'CT':-1.28, 'CC':-1.84, 'CG':-2.17,
                   'GA':-1.30, 'GT':-1.44, 'GC':-2.24, 'GG':-1.84}

        # initiation correction
        self.init_dh = 0.2
        self.init_ds = -5.7
        self.init_dg = 1.96

        # symmetry correction
        self.sym_ds = -1.4
        self.sym_dg = 0.43

        # terminal AT penalty
        self.end_dh = 2.2
        self.end_ds = 6.9
        self.end_dg = 0.05

        # 5' dangling end corrections. First base is dangling
        self.dangle5_dh = {'AA':+0.2, 'CA':+0.6, 'GA':-1.1, 'TA':-6.9,
                           'AC':-6.3, 'CC':-4.4, 'GC':-5.1, 'TC':-4.0,
                           'AG':-3.7, 'CG':-4.0, 'GG':-3.9, 'TG':-4.9,
                           'AT':-2.9, 'CT':-4.1, 'GT':-4.2, 'TT':-0.2}

        self.dangle5_dg = {'AA':-0.51, 'CA':-0.42, 'GA':-0.62, 'TA':-0.71,
                           'AC':-0.96, 'CC':-0.52, 'GC':-0.72, 'TC':-0.58,
                           'AG':-0.58, 'CG':-0.34, 'GG':-0.56, 'TG':-0.61,
                           'AT':-0.50, 'CT':-0.02, 'GT':+0.48, 'TT':-0.10}

        self.dangle5_ds = {k:((self.dangle5_dh[k] - self.dangle5_dg[k]) /
                               self.Tref) for k in self.pairs} 

        # 3' dangling end corrections. Second base is dangling
        self.dangle3_dh = {'AA':-0.5, 'AC':+4.7, 'AG':-4.1, 'AT':-3.8,
                           'CA':-5.9, 'CC':-2.6, 'CG':-3.2, 'CT':-5.2,
                           'GA':-2.1, 'GC':-0.2, 'GG':-3.9, 'GT':-4.4,
                           'TA':-0.7, 'TC':+4.4, 'TG':-1.6, 'TT':+2.9}

        self.dangle3_dg = {'AA':-0.12, 'AC':+0.28, 'AG':-0.01, 'AT':-0.13,
                           'CA':-0.82, 'CC':-0.31, 'CG':-0.01, 'CT':-0.52,
                           'GA':-0.92, 'GC':-0.23, 'GG':-0.44, 'GT':-0.35,
                           'TA':-0.48, 'TC':-0.19, 'TG':-0.50, 'TT':-0.29}

        self.dangle3_ds = {k:((self.dangle3_dh[k] - self.dangle3_dg[k]) /
                               self.Tref) for k in self.pairs} 
        
    def calc_thermo(self, duplex, T=37+273.15, saltc=None):
        """Calculate thermodynamic parameters of a DNA duplex

        Parameters
        ----------
        duplex : Duplex object
            duplex of interest
        T : float, optional
            Temperature (in K) at which to calculate Delta_G.  Default
            is 37 C (310.15 K)
        saltc : float, optional
            Concentration of monovalent salt, in mol/L (default: 1 M)

        Returns
        -------
        dH, dS, dG : tuple of float
            Delta_H, Delta_S, Delta_G of hybridization, in kcal/mol
            for dH and dG and kcal/mol/K for dS

        Notes
        -----
        Uses nearest-neighbor model to calculate enthalpy, entropy,
        and free-energy of duplex formation, using nearest-neighbor
        parameters from [1]_.  Includes corrections for symmetry,
        dangling ends, and AT termination, but does not (yet) include
        loop corrections or corrections for mismatches.  Therefore it
        implicitly assumes that there are no mismatches.

        cseq, if specified, is used only to determine dangling-end
        corrections.  The function does *not* check if it is actually
        complementary. 

        .. [1] SantaLucia and Hicks, Annu. Rev. Biophys. Biomol. Struct.
           2004.

        Examples
        --------
        In the examples below, dmodel=DNAModel()
        1) 'GTCTACC' -- a basic sequence

        >>> dH, dS, dG = dmodel.calc_thermo(Duplex('GTCTACC'))
        >>> print dH, dS
        -47.8 -0.1349

        2) 'GTCTACCA' -- a sequence with 1 end correction

        >>> dH, dS, dG = dmodel.calc_thermo(Duplex('GTCTACCA'))
        >>> print dH, dS
        -54.1 -0.1507

        3) 'GTCTAGAC' -- a self-complementary sequence

        >>> dH, dS, dG = dmodel.calc_thermo(Duplex('GTCTAGAC'))
        >>> print dH, dS
        -55.8 -0.1596
        
        4) 'TGTCTAGACA' -- a self-complementary sequence with 2 end
        corrections

        >>> dH, dS, dG = dmodel.calc_thermo(Duplex('TGTCTAGACA'))
        >>> print dH, dS
        -68.4 -0.1912

        5) 'GTCTACC' with salt correction to 0.5 M

        >>> dH, dS, dG = dmodel.calc_thermo(Duplex('GTCTACC'), saltc=0.5)
        >>> print "{:.1f} {:.4f} {:.2f}".format(dH, dS, dG)
        -47.8 -0.1364 -5.49
        
        """
        dH = dS = dG37 = 0
        s = duplex.seq
        cs = duplex.cseq
        # apply dangling end corrections first
        if str(s) != str(cs.reverse_complement()):
            print duplex
            if duplex.seq.startswith('-'):
                # add a 3' dangling end correction
                pair = str(cs[-2:])
                dH += self.dangle3_dh[pair]
                dS += self.dangle3_ds[pair]
                dG37 += self.dangle3_dg[pair]
            if cs.startswith('-'):
                # add a 3' dangling end correction
                pair = str(duplex.seq[-2:])
                dH += self.dangle3_dh[pair]
                dS += self.dangle3_ds[pair]
                dG37 += self.dangle3_dg[pair]
                # trim sequence to avoid overcounting pairs
                s = s[0:-1]
            if duplex.seq.endswith('-'):
                # add a 5' dangling end correction
                pair = str(cs[:2])
                dH += self.dangle5_dh[pair]
                dS += self.dangle5_ds[pair]
                dG37 += self.dangle5_dg[pair]
            if cs.endswith('-'):
                # add a 5' dangling end correction
                pair = str(duplex.seq[:2])
                dH += self.dangle5_dh[pair]
                dS += self.dangle5_ds[pair]
                dG37 += self.dangle5_dg[pair]
                # trim sequence to avoid overcounting pairs
                s = s[1:]
            new_duplex = Duplex(s, cs)
            print new_duplex
        
        # use dictionary comprehension to count total number of each NN pair
        nncounts = {pair:_count_overlapping(pair, str(s)) for pair
                    in self.pairs}

        dH += sum(nncounts[pair]*self.dh[pair] for pair in self.pairs)
        dS += sum(nncounts[pair]*self.ds[pair] for pair in self.pairs)
        dG37 += sum(nncounts[pair]*self.dg[pair] for pair in self.pairs)
        
        # apply symmetry correction if sequence is self-complementary
        if str(s) == str(s.reverse_complement()):
            dS += self.sym_ds 
            dG37 += self.sym_dg
    
        # apply end corrections
        if (s.startswith('A') or s.startswith('T')):
            dH += self.end_dh
            dS += self.end_ds
            dG37 += self.end_dg

        if (s.endswith('A') or s.endswith('T')):
            dH += self.end_dh
            dS += self.end_ds
            dG37 += self.end_dg

        # apply initiation correction
        dH += self.init_dh
        dS += self.init_ds
        dG37 += self.init_dg
                
        # apply salt correction
        # this implicitly assumes there is no 5' terminal phosphate, as is
        # usually the case for synthetic oligonucleotides
        if saltc:
            dS = dS + (0.368*(len(s)-1)*np.log(saltc))
            dG37 = dG37 - 0.114*(len(s)-1)*np.log(saltc)

        # apply temperature correction; assuming dCp = 0
        dG = dH - T*dS/1000

        return dH, dS/1000, dG

class RDNAModel(NNModel):
    def __init__(self):
        NNModel.__init__(self)

        # Units of energy are kcal/mole.  Units of entropy are are "e.u",
        # or cal/mol/K.  Need to divide by 1000 before multiplying by T

        # RNA-DNA duplex data.
        # From Sugimoto et al., Biochemistry 34: 11211 (1995)
        # These are in terms of the DNA NN pairs, 5' to 3' (note that
        # they are listed 3' to 5' in the Sugimoto paper
        self.dh = {'TT':-7.8, 'GT':-5.9, 'CT':-9.1, 'AT':-8.3,
                   'TG':-9.0, 'GG':-9.3, 'CG':-16.3, 'AG':-7.0,
                   'TC':-5.5, 'GC':-8.0, 'CC':-12.8, 'AC':-7.8,
                   'TA':-7.8, 'GA':-8.6, 'CA':-10.4, 'AA':-11.5}

        self.ds = {'TT':-21.9, 'GT':-12.3, 'CT':-23.5, 'AT':-23.9,
                   'TG':-26.1, 'GG':-23.2, 'CG':-47.1, 'AG':-19.7,
                   'TC':-13.5, 'GC':-17.1, 'CC':-31.9, 'AC':-21.6,
                   'TA':-23.2, 'GA':-22.9, 'CA':-28.4, 'AA':-36.4}

        self.dg = {'TT':-1.0, 'GT':-2.1, 'CT':-1.8, 'AT':-0.9,
                   'TG':-0.9, 'GG':-2.1, 'CG':-1.7, 'AG':-0.9,
                   'TC':-1.3, 'GC':-2.7, 'CC':-2.9, 'AC':-1.1,
                   'TA':-0.6, 'GA':-1.5, 'CA':-1.6, 'AA':-0.2}

        # initiation correction
        self.init_dh = 1.9
        self.init_ds = -3.9
        self.init_dg = 3.1

        
    def calc_thermo(self, duplex, T=37+273.15, saltc=None):
        """Calculate thermodynamic parameters of a DNA-RNA duplex

        Parameters
        ----------
        duplex : Duplex object
            duplex of interest
        T : float, optional
            Temperature (in K) at which to calculate Delta_G.  Default
            is 37 C (310.15 K)
        saltc : float, optional
            Concentration of monovalent salt, in mol/L (default: 1 M)

        Returns
        -------
        dH, dS, dG : tuple of float
            Delta_H, Delta_S, Delta_G of hybridization, in kcal/mol
            for dH and dG and kcal/mol/K for dS

        Notes
        -----
        Uses nearest-neighbor model to calculate enthalpy, entropy,
        and free-energy of duplex formation, using nearest-neighbor
        parameters from [1]_.  Includes corrections for symmetry,
        dangling ends, and AT termination, but does not (yet) include
        loop corrections or corrections for mismatches.  Therefore it
        implicitly assumes that there are no mismatches.

        cseq, if specified, is used only to determine dangling-end
        corrections.  The function does *not* check if it is actually
        complementary. 

        .. [1] Sugimoto et al., Biochemistry 34: 11211, 1995. 

        Examples
        --------
        In the example below, rdmodel=RDNAModel()
        1) 'AGCGTAAG' -- example from p. 11214 of Sugimoto et al. [1]_ 

        >>> dH, dS, dG = rdmodel.calc_thermo(Duplex('AGCGTAAG'))
        >>> print "{:.1f}".format(dG)
        -6.0

        2) 'AAGCGTAG' -- example from p. 11214 of Sugimoto et al. [1]_ 

        >>> dH, dS, dG = rdmodel.calc_thermo(Duplex('AAGCGTAG'))
        >>> print "{:.1f}".format(dG)
        -6.0

        """
        dH = dS = dG37 = 0
        s = duplex.seq
        
        # use dictionary comprehension to count total number of each NN pair
        nncounts = {pair:_count_overlapping(pair, str(s)) for pair
                    in self.pairs}

        dH += sum(nncounts[pair]*self.dh[pair] for pair in self.pairs)
        dS += sum(nncounts[pair]*self.ds[pair] for pair in self.pairs)
        dG37 += sum(nncounts[pair]*self.dg[pair] for pair in self.pairs)
        
        # apply initiation correction
        dH += self.init_dh
        dS += self.init_ds
        dG37 += self.init_dg

        # apply temperature correction; assuming dCp = 0
        dG = dH - T*dS/1000

        return dH, dS/1000, dG

        
# run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod(extraglobs={'dmodel': DNAModel(), 'rdmodel': RDNAModel()})
