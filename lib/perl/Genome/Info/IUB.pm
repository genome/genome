package Genome::Info::IUB;

#REVIEW fdu
#Need to fix indent problem (4-space)

use strict;
use warnings;

my %iub_as_allele_array = (
    A => ['A','A'],
    C => ['C','C'],
    G => ['G','G'],
    T => ['T','T'],
    M => ['A','C'],
    K => ['G','T'],
    Y => ['C','T'],
    R => ['A','G'],
    W => ['A','T'],
    S => ['G','C'],
    D => ['A','G','T'],
    B => ['C','G','T'],
    H => ['A','C','T'],
    V => ['A','C','G'],
    N => ['A','C','G','T'],
);

my %rna_safe_iub_as_allele_array = %iub_as_allele_array;
$rna_safe_iub_as_allele_array{'U'} = 'T';

my %iub_as_string = (
    A => 'AA',
    C => 'CC',
    G => 'GG',
    T => 'TT',
    M => 'AC',
    K => 'GT',
    Y => 'CT',
    R => 'AG',
    W => 'AT',
    S => 'CG',
    D => 'AGT',
    B => 'CGT',
    H => 'ACT',
    V => 'ACG',
    N => 'ACGT',
);

our %string_as_iub = (
    'CG' => 'S',
    'GC' => 'S',
    'AA' => 'A',
    'AT' => 'W',
    'TA' => 'W',
    'TT' => 'T',
    'ACGT' => 'N',
    'CAGT' => 'N',
    'CGAT' => 'N',
    'CGTA' => 'N',
    'AGCT' => 'N',
    'GACT' => 'N',
    'GCAT' => 'N',
    'GCTA' => 'N',
    'AGTC' => 'N',
    'GATC' => 'N',
    'GTAC' => 'N',
    'GTCA' => 'N',
    'ACTG' => 'N',
    'CATG' => 'N',
    'CTAG' => 'N',
    'CTGA' => 'N',
    'ATCG' => 'N',
    'TACG' => 'N',
    'TCAG' => 'N',
    'TCGA' => 'N',
    'ATGC' => 'N',
    'TAGC' => 'N',
    'TGAC' => 'N',
    'TGCA' => 'N',
    'GT' => 'K',
    'TG' => 'K',
    'CGT' => 'B',
    'GCT' => 'B',
    'GTC' => 'B',
    'CTG' => 'B',
    'TCG' => 'B',
    'TGC' => 'B',
    'CT' => 'Y',
    'TC' => 'Y',
    'ACG' => 'V',
    'CAG' => 'V',
    'CGA' => 'V',
    'AGC' => 'V',
    'GAC' => 'V',
    'GCA' => 'V',
    'ACT' => 'H',
    'CAT' => 'H',
    'CTA' => 'H',
    'ATC' => 'H',
    'TAC' => 'H',
    'TCA' => 'H',
    'AC' => 'M',
    'CA' => 'M',
    'AGT' => 'D',
    'GAT' => 'D',
    'GTA' => 'D',
    'ATG' => 'D',
    'TAG' => 'D',
    'TGA' => 'D',
    'CC' => 'C',
    'AG' => 'R',
    'GA' => 'R',
    'GG' => 'G',
);

my %iub_to_single_base = (
    R => 'A',
    Y => 'C',
    M => 'A',
    K => 'G',
    S => 'C',
    W => 'A',
    B => 'C',
    D => 'A',
    H => 'A',
    V => 'A',
    N => 'A',
);

sub variant_alleles_for_iub {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }
    my ($ref,$iub) = @_;

   unless(defined $ref && defined $iub) {
       return;
   }

   unless($ref =~ /[ACTGN]/i) { 
       warn "Ambiguous reference bases ($ref) not currently supported";
       return;
   }
   my @alleles = iub_to_alleles($iub);
   unless(@alleles) {
       return;
   }
   my %variants;
   foreach my $allele (@alleles) {
       if($allele ne uc($ref)) {
           $variants{$allele} = 1;
       }
   }
   return sort keys %variants;
}

sub iub_for_alleles {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }
   my @alleles = @_;
   if(@alleles != 2) {
       warn "Conversion of more than 2 alleles to IUB code is currently unsupported (".scalar @alleles." passed)";
       return; 
   }

   my %iub_for_alleles = reverse %iub_as_string;
   if(exists($iub_for_alleles{uc(join("",@alleles))})) {
       return $iub_for_alleles{uc(join("",@alleles))};
   }
   elsif(exists($iub_for_alleles{uc(join("",reverse @alleles))})) {
       return $iub_for_alleles{uc(join("", reverse @alleles))};
   }
   else {
       warn "Invalid alleles @alleles";
       return;
   }
}

sub iub_to_alleles {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }
   my ($iub) = @_;
   if(defined $iub and exists($iub_as_allele_array{uc $iub})) { 
       return @{$iub_as_allele_array{uc $iub}};
   }
   else {
       return;
   }
}

sub iub_to_bases {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }
   my ($iub) = @_;
   return _iub_to_bases(\%iub_as_allele_array, $iub);
}

sub rna_safe_iub_to_bases {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }
   my ($iub) = @_;
   return _iub_to_bases(\%rna_safe_iub_as_allele_array, $iub);
}

sub _iub_to_bases {
   my %iub_as_array = %{shift(@_)};
   my ($iub) = @_;
   my %bases = map {$_ => 1} @{$iub_as_array{uc $iub}};
   return sort keys %bases;
}

sub iub_to_string {
    my $base = pop;
    $base = uc $base if $base;
    
    return $iub_as_string{$base} if $base and $iub_as_string{$base};
    return;
}


sub string_to_iub {
    my $string = pop;
    $string    = uc $string;
    my $iub = $string_as_iub{$string};
    unless ($iub) {
        Carp::confess("Cannot translate string $string into IUB code, expected: " . join(",",sort keys %string_as_iub));
    }
    return $iub;
}


sub reference_iub_to_base {
    my $class;
    if ((defined($_[0]))&&($_[0] eq __PACKAGE__)) {
        $class = shift;
    }

    my($reference_with_iub) = @_;
    $reference_with_iub = uc $reference_with_iub if $reference_with_iub;

    return $reference_with_iub unless $reference_with_iub;
    
    my $reference_without_iubs; 
    
    for my $base(split(//, $reference_with_iub)){
        if($iub_to_single_base{$base}){
            $reference_without_iubs .= $iub_to_single_base{$base};
        }else{
            $reference_without_iubs .= $base;
        }
    }
    return $reference_without_iubs;
}

1;

=pod

=head1 Name

Genome::Info::IUB

=head1 Synopsis

Genome::Info::IUB contains methods for working with IUB codes

=head1 Usage

my $iub = Genome::Info::IUB::iub_for_alleles("A","T"); #returns 'W'

my @alleles = Genome::Info::IUB::iub_to_alleles("A"); #returns ('A','A')

my @bases = Genome::Info::IUB::iub_to_bases("A"); #returns ('A')

my @variants = Genome::Info::IUB::variant_alleles_for_iub("T","W"); #returns ('A') 

=head1 Methods

=head2 iub_for_alleles 

=over 

=item I<Synopsis>

converts an allele string into an IUB code

=item I<Arguments>

a string of two bases 

=item I<Returns>

the IUB code for the two bases passed in or undef if no such IUB code exists or the alleles were invalid

=back

=head2 iub_to_alleles 

=over

=item I<Synopsis>

converts an IUB code to the alleles it represents

=item I<Arguments>

an IUB code

=item I<Returns>

a list of the alleles that the IUB code represents with at least two bases returned or undef if invalid. 

=back

=head2 iub_to_bases 

=over

=item I<Synopsis>

converts an IUB code to the bases it represents

=item I<Arguments>

an IUB code

=item I<Returns>

a list of the bases that the IUB code represent or undef if invalid. Non-ambiguous bases return that base only.

=back

=head2 variant_alleles_for_iub 

=over

=item I<Synopsis>

converts an IUB code to variants alleles given a reference base

=item I<Arguments>

a reference base and the IUB code

=item I<Returns>

a list of the variant alleles for the IUB code given the reference base or undef if invalid

=back

=head1 Author(s)

B<David Larson> I<dlarson@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Info/CodonToAminoAcid.pm $
#$Id: CodonToAminoAcid.pm 34977 2008-05-23 22:34:14Z ebelter $
