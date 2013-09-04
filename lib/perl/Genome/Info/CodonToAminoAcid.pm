package Genome::Info::CodonToAminoAcid;

use strict;
use warnings;

# From Chapter 8 codon2aa
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

my @indels = qw(
    -TA  -TC  -TG  -TT  -CA  -CC  -CG  -CT  -AC  -AT  -AA  -AG  -GA  -GC  -GG  -GT
    T-A  T-C  T-G  T-T  C-A  C-C  C-G  C-T  A-C  A-T  A-A  A-G  G-A  G-C  G-G  G-T
    TA-  TC-  TG-  TT-  CA-  CC-  CG-  CT-  AC-  AT-  AA-  AG-  GA-  GC-  GG-  GT-
);


my @refseq_allele = qw(
    +TA  +TC  +TG  +TT  +CA  +CC  +CG  +CT  +AC  +AT  +AA  +AG  +GA  +GC  +GG  +GT
    T+A  T+C  T+G  T+T  C+A  C+C  C+G  C+T  A+C  A+T  A+A  A+G  G+A  G+C  G+G  G+T
    TA+  TC+  TG+  TT+  CA+  CC+  CG+  CT+  AC+  AT+  AA+  AG+  GA+  GC+  GG+  GT+
);


my %names = (
    TCA => 'S',  # Serine
    TCC => 'S',  # Serine
    TCG => 'S',  # Serine
    TCT => 'S',  # Serine
    TTC => 'F',  # Phenylalanine
    TTT => 'F',  # Phenylalanine
    TTA => 'L',  # Leucine
    TTG => 'L',  # Leucine
    TAC => 'Y',  # Tyrosine
    TAT => 'Y',  # Tyrosine
    TGC => 'C',  # Cysteine
    TGT => 'C',  # Cysteine
    TGG => 'W',  # Tryptophan
    CTA => 'L',  # Leucine
    CTC => 'L',  # Leucine
    CTG => 'L',  # Leucine
    CTT => 'L',  # Leucine
    CCA => 'P',  # Proline
    CCC => 'P',  # Proline
    CCG => 'P',  # Proline
    CCT => 'P',  # Proline
    CAC => 'H',  # Histidine
    CAT => 'H',  # Histidine
    CAA => 'Q',  # Glutamine
    CAG => 'Q',  # Glutamine
    CGA => 'R',  # Arginine
    CGC => 'R',  # Arginine
    CGG => 'R',  # Arginine
    CGT => 'R',  # Arginine
    ATA => 'I',  # Isoleucine
    ATC => 'I',  # Isoleucine
    ATT => 'I',  # Isoleucine
    ATG => 'M',  # Methionine
    ACA => 'T',  # Threonine
    ACC => 'T',  # Threonine
    ACG => 'T',  # Threonine
    ACT => 'T',  # Threonine
    AAC => 'N',  # Asparagine
    AAT => 'N',  # Asparagine
    AAA => 'K',  # Lysine
    AAG => 'K',  # Lysine
    AGC => 'S',  # Serine
    AGT => 'S',  # Serine
    AGA => 'R',  # Arginine
    AGG => 'R',  # Arginine
    GTA => 'V',  # Valine
    GTC => 'V',  # Valine
    GTG => 'V',  # Valine
    GTT => 'V',  # Valine
    GCA => 'A',  # Alanine
    GCC => 'A',  # Alanine
    GCG => 'A',  # Alanine
    GCT => 'A',  # Alanine
    GAC => 'D',  # Aspartic Acid
    GAT => 'D',  # Aspartic Acid
    GAA => 'E',  # Glutamic Acid
    GAG => 'E',  # Glutamic Acid
    GGA => 'G',  # Glycine
    GGC => 'G',  # Glycine
    GGG => 'G',  # Glycine
    GGT => 'G',  # Glycine
    XTA => 'Z',  #Discrepant Genotypes in Overlapping Data
    XTC => 'Z',  #Discrepant Genotypes in Overlapping Data
    XTG => 'Z',  #Discrepant Genotypes in Overlapping Data
    XTT => 'Z',  #Discrepant Genotypes in Overlapping Data
    XCA => 'Z',  #Discrepant Genotypes in Overlapping Data
    XCC => 'Z',  #Discrepant Genotypes in Overlapping Data
    XCG => 'Z',  #Discrepant Genotypes in Overlapping Data
    XCT => 'Z',  #Discrepant Genotypes in Overlapping Data
    XAC => 'Z',  #Discrepant Genotypes in Overlapping Data
    XAT => 'Z',  #Discrepant Genotypes in Overlapping Data
    XAA => 'Z',  #Discrepant Genotypes in Overlapping Data
    XAG => 'Z',  #Discrepant Genotypes in Overlapping Data
    XGA => 'Z',  #Discrepant Genotypes in Overlapping Data
    XGC => 'Z',  #Discrepant Genotypes in Overlapping Data
    XGG => 'Z',  #Discrepant Genotypes in Overlapping Data
    XGT => 'Z',  #Discrepant Genotypes in Overlapping Data
    TXA => 'Z',  #Discrepant Genotypes in Overlapping Data
    TXC => 'Z',  #Discrepant Genotypes in Overlapping Data
    TXG => 'Z',  #Discrepant Genotypes in Overlapping Data
    TXT => 'Z',  #Discrepant Genotypes in Overlapping Data
    CXA => 'Z',  #Discrepant Genotypes in Overlapping Data
    CXC => 'Z',  #Discrepant Genotypes in Overlapping Data
    CXG => 'Z',  #Discrepant Genotypes in Overlapping Data
    CXT => 'Z',  #Discrepant Genotypes in Overlapping Data
    AXC => 'Z',  #Discrepant Genotypes in Overlapping Data
    AXT => 'Z',  #Discrepant Genotypes in Overlapping Data
    AXA => 'Z',  #Discrepant Genotypes in Overlapping Data
    AXG => 'Z',  #Discrepant Genotypes in Overlapping Data
    GXA => 'Z',  #Discrepant Genotypes in Overlapping Data
    GXC => 'Z',  #Discrepant Genotypes in Overlapping Data
    GXG => 'Z',  #Discrepant Genotypes in Overlapping Data
    GXT => 'Z',  #Discrepant Genotypes in Overlapping Data
    TAX => 'Z',  #Discrepant Genotypes in Overlapping Data
    TCX => 'S',  #Discrepant Genotypes in Overlapping Data #Z->S Ser
    TGX => 'Z',  #Discrepant Genotypes in Overlapping Data
    TTX => 'Z',  #Discrepant Genotypes in Overlapping Data
    CAX => 'Z',  #Discrepant Genotypes in Overlapping Data
    CCX => 'P',  #Discrepant Genotypes in Overlapping Data #Z->P Pro
    CGX => 'R',  #Discrepant Genotypes in Overlapping Data #Z->R Arg
    CTX => 'L',  #Discrepant Genotypes in Overlapping Data #Z->L Leu
    ACX => 'T',  #Discrepant Genotypes in Overlapping Data #Z->T Thr
    ATX => 'Z',  #Discrepant Genotypes in Overlapping Data
    AAX => 'Z',  #Discrepant Genotypes in Overlapping Data
    AGX => 'Z',  #Discrepant Genotypes in Overlapping Data
    GAX => 'Z',  #Discrepant Genotypes in Overlapping Data
    GCX => 'A',  #Discrepant Genotypes in Overlapping Data #Z->A Ala
    GGX => 'G',  #Discrepant Genotypes in Overlapping Data #Z->G Gly
    GTX => 'V',  #Discrepant Genotypes in Overlapping Data #Z->V Val
);

#stop codons
my %stop_codons = (
    TAA => ['X', 'OCH'],
    TAG => ['X', 'AMB'],
    TGA => ['X', 'OPA'],
);

my %convert = (
    S => 'Ser',
    F => 'Phe',
    L => 'Leu',
    Y => 'Tyr',
    C => 'Cys',
    W => 'Trp',
    P => 'Pro',
    H => 'His',
    Q => 'Gln',
    R => 'Arg',
    I => 'Ile',
    M => 'Met',
    T => 'Thr',
    N => 'Asn',
    K => 'Lys',
    V => 'Val',
    A => 'Ala',
    D => 'Asp',
    E => 'Glu',
    G => 'Gly',
    Z => 'Z',
);


sub single_letter {
    my %_stop_codons = _stop_codons('single');
    my %_others      = _others();
    return (%names, %_stop_codons, %_others);
}


sub three_letter {
    my %_stop_codons = _stop_codons('three');
    my %_others      = _others();
    my %three_letter;
    
    for my $codon (keys %names) {
        my $single = $names{$codon};
        my $three  = $convert{$single};
        unless ($three) {
            die "There is no three letter amino acid name for codon: $codon single letter: $single\n";
        }
        $three_letter{$codon} = $three;
    }
    return (%three_letter, %_stop_codons, %_others);
}

sub _others {
    my %indels;
    map{$indels{$_}= 'indel'}@indels;

    my %refseq_allele;
    map{$refseq_allele{$_}='refseq allele'}@refseq_allele;

    return (%indels, %refseq_allele);
}

sub _stop_codons {
    my $type = pop;
    my $pos  = $type eq 'single' ? 0 : 1;

    my %_stop_codons;
    for my $codon (keys %stop_codons) {
        $_stop_codons{$codon} = $stop_codons{$codon}->[$pos];
    }

    return %_stop_codons;
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
