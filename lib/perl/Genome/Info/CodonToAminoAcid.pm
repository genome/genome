package Genome::Info::CodonToAminoAcid;

#REVIEW fdu
#1. No need to pull hash table outside. Just return the hash inside the subroutine.
#2. Couldn't find any other module using this under Genome tree

use strict;
use warnings;

# From Chapter 8 codon2aa
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup
my %single_letter = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'X',    # Stop
    'TAG' => 'X',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'X',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    '-TA' => 'indel', #Indel
    '-TC' => 'indel', #Indel
    '-TG' => 'indel', #Indel
    '-TT' => 'indel', #Indel
    '-CA' => 'indel', #Indel
    '-CC' => 'indel', #Indel
    '-CG' => 'indel', #Indel
    '-CT' => 'indel', #Indel
    '-AC' => 'indel', #Indel
    '-AT' => 'indel', #Indel
    '-AA' => 'indel', #Indel
    '-AG' => 'indel', #Indel
    '-GA' => 'indel', #Indel
    '-GC' => 'indel', #Indel
    '-GG' => 'indel', #Indel
    '-GT' => 'indel', #Indel
    'T-A' => 'indel', #Indel
    'T-C' => 'indel', #Indel
    'T-G' => 'indel', #Indel
    'T-T' => 'indel', #Indel
    'C-A' => 'indel', #Indel
    'C-C' => 'indel', #Indel
    'C-G' => 'indel', #Indel
    'C-T' => 'indel', #Indel
    'A-C' => 'indel', #Indel
    'A-T' => 'indel', #Indel
    'A-A' => 'indel', #Indel
    'A-G' => 'indel', #Indel
    'G-A' => 'indel', #Indel
    'G-C' => 'indel', #Indel
    'G-G' => 'indel', #Indel
    'G-T' => 'indel', #Indel
    'TA-' => 'indel', #Indel
    'TC-' => 'indel', #Indel
    'TG-' => 'indel', #Indel
    'TT-' => 'indel', #Indel
    'CA-' => 'indel', #Indel
    'CC-' => 'indel', #Indel
    'CG-' => 'indel', #Indel
    'CT-' => 'indel', #Indel
    'AC-' => 'indel', #Indel
    'AT-' => 'indel', #Indel
    'AA-' => 'indel', #Indel
    'AG-' => 'indel', #Indel
    'GA-' => 'indel', #Indel
    'GC-' => 'indel', #Indel
    'GG-' => 'indel', #Indel
    'GT-' => 'indel', #Indel
    '+TA' => 'refseq allele', #No Indel
    '+TC' => 'refseq allele', #No Indel
    '+TG' => 'refseq allele', #No Indel
    '+TT' => 'refseq allele', #No Indel
    '+CA' => 'refseq allele', #No Indel
    '+CC' => 'refseq allele', #No Indel
    '+CG' => 'refseq allele', #No Indel
    '+CT' => 'refseq allele', #No Indel
    '+AC' => 'refseq allele', #No Indel
    '+AT' => 'refseq allele', #No Indel
    '+AA' => 'refseq allele', #No Indel
    '+AG' => 'refseq allele', #No Indel
    '+GA' => 'refseq allele', #No Indel
    '+GC' => 'refseq allele', #No Indel
    '+GG' => 'refseq allele', #No Indel
    '+GT' => 'refseq allele', #No Indel
    'T+A' => 'refseq allele', #No Indel
    'T+C' => 'refseq allele', #No Indel
    'T+G' => 'refseq allele', #No Indel
    'T+T' => 'refseq allele', #No Indel
    'C+A' => 'refseq allele', #No Indel
    'C+C' => 'refseq allele', #No Indel
    'C+G' => 'refseq allele', #No Indel
    'C+T' => 'refseq allele', #No Indel
    'A+C' => 'refseq allele', #No Indel
    'A+T' => 'refseq allele', #No Indel
    'A+A' => 'refseq allele', #No Indel
    'A+G' => 'refseq allele', #No Indel
    'G+A' => 'refseq allele', #No Indel
    'G+C' => 'refseq allele', #No Indel
    'G+G' => 'refseq allele', #No Indel
    'G+T' => 'refseq allele', #No Indel
    'TA+' => 'refseq allele', #No Indel
    'TC+' => 'refseq allele', #No Indel
    'TG+' => 'refseq allele', #No Indel
    'TT+' => 'refseq allele', #No Indel
    'CA+' => 'refseq allele', #No Indel
    'CC+' => 'refseq allele', #No Indel
    'CG+' => 'refseq allele', #No Indel
    'CT+' => 'refseq allele', #No Indel
    'AC+' => 'refseq allele', #No Indel
    'AT+' => 'refseq allele', #No Indel
    'AA+' => 'refseq allele', #No Indel
    'AG+' => 'refseq allele', #No Indel
    'GA+' => 'refseq allele', #No Indel
    'GC+' => 'refseq allele', #No Indel
    'GG+' => 'refseq allele', #No Indel
    'GT+' => 'refseq allele', #No Indel
    'XTA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TCX' => 'S', #Discrepant Genotypes in Overlapping Data #Z->S Ser
    'TGX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TTX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CCX' => 'P', #Discrepant Genotypes in Overlapping Data #Z->P Pro
    'CGX' => 'R', #Discrepant Genotypes in Overlapping Data #Z->R Arg
    'CTX' => 'L', #Discrepant Genotypes in Overlapping Data #Z->L Leu
    'ACX' => 'T', #Discrepant Genotypes in Overlapping Data #Z->T Thr
    'ATX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AGX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GCX' => 'A', #Discrepant Genotypes in Overlapping Data #Z->A Ala
    'GGX' => 'G', #Discrepant Genotypes in Overlapping Data #Z->G Gly
    'GTX' => 'V', #Discrepant Genotypes in Overlapping Data #Z->V Val
    );

sub single_letter {
    return %single_letter;
}

# From Chapter 8 codon2aa
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup
my %three_letter = (
    'TCA' => 'Ser',    # Serine
    'TCC' => 'Ser',    # Serine
    'TCG' => 'Ser',    # Serine
    'TCT' => 'Ser',    # Serine
    'TTC' => 'Phe',    # Phenylalanine
    'TTT' => 'Phe',    # Phenylalanine
    'TTA' => 'Leu',    # Leucine
    'TTG' => 'Leu',    # Leucine
    'TAC' => 'Tyr',    # Tyrosine
    'TAT' => 'Tyr',    # Tyrosine
    'TAA' => 'OCH',    # Stop
    'TAG' => 'AMB',    # Stop
    'TGC' => 'Cys',    # Cysteine
    'TGT' => 'Cys',    # Cysteine
    'TGA' => 'OPA',    # Stop
    'TGG' => 'Trp',    # Tryptophan
    'CTA' => 'Leu',    # Leucine
    'CTC' => 'Leu',    # Leucine
    'CTG' => 'Leu',    # Leucine
    'CTT' => 'Leu',    # Leucine
    'CCA' => 'Pro',    # Proline
    'CCC' => 'Pro',    # Proline
    'CCG' => 'Pro',    # Proline
    'CCT' => 'Pro',    # Proline
    'CAC' => 'His',    # Histidine
    'CAT' => 'His',    # Histidine
    'CAA' => 'Gln',    # Glutamine
    'CAG' => 'Gln',    # Glutamine
    'CGA' => 'Arg',    # Arginine
    'CGC' => 'Arg',    # Arginine
    'CGG' => 'Arg',    # Arginine
    'CGT' => 'Art',    # Arginine
    'ATA' => 'Ile',    # Isoleucine
    'ATC' => 'Ile',    # Isoleucine
    'ATT' => 'Ile',    # Isoleucine
    'ATG' => 'Met',    # Methionine
    'ACA' => 'Thr',    # Threonine
    'ACC' => 'Thr',    # Threonine
    'ACG' => 'Thr',    # Threonine
    'ACT' => 'Thr',    # Threonine
    'AAC' => 'Asn',    # Asparagine
    'AAT' => 'Asn',    # Asparagine
    'AAA' => 'Lys',    # Lysine
    'AAG' => 'Lys',    # Lysine
    'AGC' => 'Ser',    # Serine
    'AGT' => 'Ser',    # Serine
    'AGA' => 'Arg',    # Arginine
    'AGG' => 'Arg',    # Arginine
    'GTA' => 'Val',    # Valine
    'GTC' => 'Val',    # Valine
    'GTG' => 'Val',    # Valine
    'GTT' => 'Val',    # Valine
    'GCA' => 'Ala',    # Alanine
    'GCC' => 'Ala',    # Alanine
    'GCG' => 'Ala',    # Alanine
    'GCT' => 'Ala',    # Alanine
    'GAC' => 'Asp',    # Aspartic Acid
    'GAT' => 'Asp',    # Aspartic Acid
    'GAA' => 'Glu',    # Glutamic Acid
    'GAG' => 'Glu',    # Glutamic Acid
    'GGA' => 'Gly',    # Glycine
    'GGC' => 'Gly',    # Glycine
    'GGG' => 'Gly',    # Glycine
    'GGT' => 'Gly',    # Glycine
    '-TA' => 'indel', #Indel
    '-TC' => 'indel', #Indel
    '-TG' => 'indel', #Indel
    '-TT' => 'indel', #Indel
    '-CA' => 'indel', #Indel
    '-CC' => 'indel', #Indel
    '-CG' => 'indel', #Indel
    '-CT' => 'indel', #Indel
    '-AC' => 'indel', #Indel
    '-AT' => 'indel', #Indel
    '-AA' => 'indel', #Indel
    '-AG' => 'indel', #Indel
    '-GA' => 'indel', #Indel
    '-GC' => 'indel', #Indel
    '-GG' => 'indel', #Indel
    '-GT' => 'indel', #Indel
    'T-A' => 'indel', #Indel
    'T-C' => 'indel', #Indel
    'T-G' => 'indel', #Indel
    'T-T' => 'indel', #Indel
    'C-A' => 'indel', #Indel
    'C-C' => 'indel', #Indel
    'C-G' => 'indel', #Indel
    'C-T' => 'indel', #Indel
    'A-C' => 'indel', #Indel
    'A-T' => 'indel', #Indel
    'A-A' => 'indel', #Indel
    'A-G' => 'indel', #Indel
    'G-A' => 'indel', #Indel
    'G-C' => 'indel', #Indel
    'G-G' => 'indel', #Indel
    'G-T' => 'indel', #Indel
    'TA-' => 'indel', #Indel
    'TC-' => 'indel', #Indel
    'TG-' => 'indel', #Indel
    'TT-' => 'indel', #Indel
    'CA-' => 'indel', #Indel
    'CC-' => 'indel', #Indel
    'CG-' => 'indel', #Indel
    'CT-' => 'indel', #Indel
    'AC-' => 'indel', #Indel
    'AT-' => 'indel', #Indel
    'AA-' => 'indel', #Indel
    'AG-' => 'indel', #Indel
    'GA-' => 'indel', #Indel
    'GC-' => 'indel', #Indel
    'GG-' => 'indel', #Indel
    'GT-' => 'indel', #Indel
    '+TA' => 'refseq allele', #No Indel
    '+TC' => 'refseq allele', #No Indel
    '+TG' => 'refseq allele', #No Indel
    '+TT' => 'refseq allele', #No Indel
    '+CA' => 'refseq allele', #No Indel
    '+CC' => 'refseq allele', #No Indel
    '+CG' => 'refseq allele', #No Indel
    '+CT' => 'refseq allele', #No Indel
    '+AC' => 'refseq allele', #No Indel
    '+AT' => 'refseq allele', #No Indel
    '+AA' => 'refseq allele', #No Indel
    '+AG' => 'refseq allele', #No Indel
    '+GA' => 'refseq allele', #No Indel
    '+GC' => 'refseq allele', #No Indel
    '+GG' => 'refseq allele', #No Indel
    '+GT' => 'refseq allele', #No Indel
    'T+A' => 'refseq allele', #No Indel
    'T+C' => 'refseq allele', #No Indel
    'T+G' => 'refseq allele', #No Indel
    'T+T' => 'refseq allele', #No Indel
    'C+A' => 'refseq allele', #No Indel
    'C+C' => 'refseq allele', #No Indel
    'C+G' => 'refseq allele', #No Indel
    'C+T' => 'refseq allele', #No Indel
    'A+C' => 'refseq allele', #No Indel
    'A+T' => 'refseq allele', #No Indel
    'A+A' => 'refseq allele', #No Indel
    'A+G' => 'refseq allele', #No Indel
    'G+A' => 'refseq allele', #No Indel
    'G+C' => 'refseq allele', #No Indel
    'G+G' => 'refseq allele', #No Indel
    'G+T' => 'refseq allele', #No Indel
    'TA+' => 'refseq allele', #No Indel
    'TC+' => 'refseq allele', #No Indel
    'TG+' => 'refseq allele', #No Indel
    'TT+' => 'refseq allele', #No Indel
    'CA+' => 'refseq allele', #No Indel
    'CC+' => 'refseq allele', #No Indel
    'CG+' => 'refseq allele', #No Indel
    'CT+' => 'refseq allele', #No Indel
    'AC+' => 'refseq allele', #No Indel
    'AT+' => 'refseq allele', #No Indel
    'AA+' => 'refseq allele', #No Indel
    'AG+' => 'refseq allele', #No Indel
    'GA+' => 'refseq allele', #No Indel
    'GC+' => 'refseq allele', #No Indel
    'GG+' => 'refseq allele', #No Indel
    'GT+' => 'refseq allele', #No Indel
    'XTA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XTT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XCT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XAG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'XGT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXA' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXC' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXG' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GXT' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TCX' => 'Ser', #Discrepant Genotypes in Overlapping Data
    'TGX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'TTX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'CCX' => 'Pro', #Discrepant Genotypes in Overlapping Data
    'CGX' => 'Arg', #Discrepant Genotypes in Overlapping Data
    'CTX' => 'Lue', #Discrepant Genotypes in Overlapping Data
    'ACX' => 'Ale', #Discrepant Genotypes in Overlapping Data
    'ATX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'AGX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GAX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GCX' => 'Z', #Discrepant Genotypes in Overlapping Data
    'GGX' => 'Gly', #Discrepant Genotypes in Overlapping Data
    'GTX' => 'Val', #Discrepant Genotypes in Overlapping Data
    );

sub three_letter {
    return %three_letter;
}

1;

=pod

############################################################
#
# GSC_AminoAcidTranslation.pm borrowed from BeginPerlBioinfo.pm
# - a library of subroutines
#   from the examples and text in the book:
#
# Beginning Perl for Bioinformatics
# by James Tisdall
#
# published by O'Reilly & Associates
# (c) 2001 James Tisdall
#
# to your code, making sure the module is in the same directory
# as the program you are using it from, or another place where
# Perl can find it (see the discussion in the book for other locations).
#
# Version 20011230
#   incorporates a few errata and bug fixes
# 
# mmclella for polyphred type indel genotypes including (+/-/X)
#          where + indicates a match to refseq allle
#          where - indicates an insertion or a deletion
#          where X indicates discrepant/conflicting data.
#
############################################################

=cut

#$HeadURL$
#$Id$
