# review gsanders
# This and HugoGene.pm could probably use some work...
# These are currently only used for Genome:Model:Tools:Somatic::UcscAnnotator.pm.
# Previously there were just home-directory modules. We could trim these down to just the necessary code for ucsc...

package Genome::Utility::HugoGene::HugoGene;

use strict;
use Carp;
use Genome;

# Used to hold data from Hugo web site

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};

    # Unique values
    $self->{ID} = undef;             # Unique ID number (HGNC ID after removing 'HGNC:')
    $self->{SYMBOL} = undef;         # HUGO symbol (gene name)
    $self->{NAME} = undef;           # Gene description (string with white space)
    $self->{TYPE} = undef;           # Locus type ('pseudogene', 'gene with protein product, inferred', etc.)
    $self->{CHR} = undef;            # Chromosome assignment (includes band)
    $self->{UNIPROT} = undef;        # Uniprot name (mapped data supplied by UniProt)
    $self->{ENSEMBL} = undef;        # Ensembl name (mapped data supplied by Ensembl)
    $self->{OMIM} = undef;           # Online Mendelian Inheritance in Man name (mapped data supplied by NCBI)
    $self->{UCSC} = undef;           # UCSC ID (mapped data supplied by UCSC)
    $self->{START} = undef;          # Nucleotide start on chr (not supplied by Hugo)
    $self->{STOP} = undef;           # Nucleotide end on chr (not supplied by Hugo) 

    # Not necessarily unique
    $self->{ENTREZ} = ();            # Hash of Entrez (LocusId) name(s) 
    $self->{REFSEQ} = ();            # Hash of Ref Seq name(s) 
    $self->{VEGA} = ();              # Hash of Vega name(s) 
    $self->{PREVIOUS_SYMBOLS} = ();  # Hash of previous Hugo symbols 
                                     # (previous symbol may associate with more than one unique Hugo ID/symbol)
    $self->{ALIAS} = ();             # Hash of aliases 
                                     # (an alias may associate with more than one unique Hugo ID/symbol) 

    bless ($self, $class);
    return $self;
}


### Set/Get routines for all unique variables
sub id {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{ID} = shift; $self->{ID} =~ s/\s+//g; }
    return $self->{ID};
}

sub symbol {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{SYMBOL} = shift; $self->{SYMBOL} =~ s/\s+//g;   }
    return $self->{SYMBOL};
}
    
sub name {
    # Whitespace expected
    my $self = shift;
    if (@_) { $self->{NAME} = shift; }
    return $self->{NAME};
}
  
sub type {
    # Whitespace expected
    my $self = shift;
    if (@_) { $self->{TYPE} = shift; }
    return $self->{TYPE};
}

sub chr {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
	$self->{CHR} = shift; 
	$self->{CHR} =~ s/\s+//g;  
    }
    return $self->{CHR};
}

sub uniprot { 
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{UNIPROT} = shift; $self->{UNIPROT} =~ s/\s+//g;  }
    return $self->{UNIPROT};
}
    
sub ensembl {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{ENSEMBL} = shift; $self->{ENSEMBL} =~ s/\s+//g;  }
    return $self->{ENSEMBL};
}
      
sub omim {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{OMIM} = shift; $self->{OMIM} =~ s/\s+//g;  }
    return $self->{OMIM};
}   

sub ucsc {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{UCSC} = shift; $self->{UCSC} =~ s/\s+//g;  }
    return $self->{UCSC};
}

sub start {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{START} = shift; $self->{START} =~ s/\s+//g;  }
    return $self->{START};
}
     
sub stop {
    # Don't want whitespace
    my $self = shift;
    if (@_) { $self->{STOP} = shift; $self->{STOP} =~ s/\s+//g;  }
    return $self->{STOP};
}
      
    

### Set routines for variables stored in hash

sub addEntrezName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{ENTREZ}}{$id} = 1;
}

sub addRefSeqName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{REFSEQ}}{$id} = 1;
}

sub addVegaName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{VEGA}}{$id} = 1;
}


sub addPreviousSymbol {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{PREVIOUS_SYMBOLS}}{$id} = 1;
}

sub addAlias {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{ALIAS}}{$id} = 1;
}
    

## Get routines for variables stored in hash
# All return ref to hash 

sub getAllEntrezNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{ENTREZ}};
}
 
sub getAllRefSeqNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{REFSEQ}};
}
 
sub getAllVegaNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{VEGA}};
}

sub getAllPreviousSymbols {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{PREVIOUS_SYMBOLS}};
}


sub getAllAliases {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{ALIAS}};
}



## Methods to see if this Hugo symbol refers to the same
## gene as identified 
sub isPreviousSymbol {
    # Input: Hugo name
    # Return: 1 if the input name matches one of the
    #         previous symbols of this object

    my ($self, $hugoName) = @_;
    $hugoName =~ s/\s+//g;
    foreach ( keys %{$self->getAllPreviousSymbols()} ) {
	if ( $_ eq $hugoName ) { return 1; }
    }

    return 0;
}


sub isAlias {
    # Input: Hugo name
    # Return: 1 if the input name matches one of the
    #         alias symbols of this object

    my ($self, $hugoName) = @_;
    $hugoName =~ s/\s+//g;
    foreach ( keys %{$self->getAllAliases()} ) {
	if ( $_ eq $hugoName ) { return 1; }
    }

    return 0;
}

sub hasVega {
    # Input: Vega gene ID
    # Return: 1 if this object has the Vega gene ID

    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllVegaNames()} ) {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;
}

sub hasEntrez {
    # Input: Entrez gene ID
    # Return: 1 if this object has the Entrez gene ID

    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllEntrezNames()} ) {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;

}

sub hasRefSeq {
    # Input: RefSeq gene ID
    # Return: 1 if this object has the RefSeq gene ID

    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllRefSeqNames()} ) {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;

}

sub overlaps {
    # Input: start, stop, tolerance (optional)
    # Return: 1 if the start, stop overlap the start, stop of this gene
    #         within 'tolerance' base pairs
    # Returns undef if this gene does not have start, stop
    # ASSUMES: we are talking about the same chromosome.  Does not check 
    # chromosome.  

    my ( $self, $start, $stop, $tolerance ) = @_;
    if ( !defined $tolerance ) { $tolerance = 0; }
    if ( !defined $self->{START} || !defined $self->{STOP} ) {
	carp "WARNING: $self->{SYMBOL} does not have start, stop.  Can not determine overlap \n";
	return undef;
    }

    # Make sure start < stop for both
    my $thisStart = $self->start();
    my $thisStop = $self->stop();
    if ( $thisStart > $thisStop ) { ($thisStart, $thisStop) = ($thisStop, $thisStart); }
    if ( $start > $stop ) { ($start, $stop) = ($stop, $start); }

    return ( $stop + $tolerance > $thisStart && $start - $tolerance < $thisStop );
}

    



sub sameMajorChrBand {
    # Input: chromosome (in any format)
    # Return: 1 if this object is on the same major chromosome band
    # i.e. It does not look at the number after the decimal point
    # 7p14.20 ~ 7p14 ~ 7p14.2  etc.
    # If one band is 'ter' or 'cen', it will not necessarily return the
    # correct answer since the 'ter' band can have a standard name

    my ( $self, $inputChr ) = @_;
    my ( $chr, $arm, $band, $inputArm, $inputBand );
    $chr = $self->chr();

    # This should not be called if $chr is not defined
    (defined $chr) ||
	carp $self->symbol(), " does not have a defined chromosome";
  
    # Remove common junk
    $chr =~ s/chr//i; $inputChr =~ s/chr//i;
    $chr =~ s/\_//;  $inputChr =~ s/\_//;
    $chr =~ s/\s+//g; $inputChr =~ s/\s+//g;

    # Change to uppercase in case chromosomes are X or Y
    $chr = uc $chr;
    $inputChr = uc $inputChr;

    # Return 0 if the chromosomes are different
    my ($entireChr, $entireInputChr );
    if ( $chr =~ /(X|Y|\d+)(p|q|cen|ter)?/i ) { $entireChr = $1; }
    if ( $inputChr =~ /(X|Y|\d+)(p|q|cen|ter)?/i ) { $entireInputChr = $1; }
    if ( !defined $entireChr || !defined $entireInputChr ) {
	carp "Unexpected format for \$inputChr: '$inputChr' sent to $self->{SYMBOL} or \$chr unexpected format";
	return 0;
    }	    
    if ( $entireChr ne $entireInputChr ) { return 0; }

    # One or both bands might be expressed as a range.
    

#### maybe make a helper function to look at two bands......
# so it can be called with potentially 4 different combinations....






    # They might match now
    if ( $chr eq $inputChr ) { return 1; }

    # If one of the chromosomes has cen, pter, or qter,
    # just go by chromosome number
    if ( $chr =~ /ter/i || $inputChr =~ /ter/i ||
	 $chr =~ /cen/i || $inputChr =~ /cen/i ) {

	if ( $chr =~ /(\w+)(p|q|cen)/i ) { $chr = $1; }
	if ( $inputChr =~ /(\w+)(p|q|cen)/i ) { $inputChr = $1; }
	return ( $inputChr eq $chr );
    }
    
    # Chromosomes can be in format '$chr[pq]\d+.\d+'
    # or can be given as a range 
    if ( $chr =~ /(\w+)([pq])(\d*)?/i ) {
	$chr = $1; $arm = $2; $band = $3;
    }
    if ( $inputChr =~ /(\w+)([pq])(\d*)?/i ) {
	$inputChr = $1; $inputArm = $2; $inputBand = $3;
    }

       
    # Return 0 if chromosomes are different
    if ( $inputChr ne $chr ) { return 0; }

    # If one of the chromosomes does not have an arm,
    # it is OK to return 1 since the chromosomes are the same
    if ( !defined $arm || !defined $inputArm ||
	 $arm eq "" || $inputArm eq "" ) { return 1; }
    
    # If the arms are different, return 0
    if ( $arm ne $inputArm ) { return 0; }
    
    # The chromosomes and arms are the same, allow the
    # bands to differ by 1
    if ( defined $band && defined $inputBand && $band ne "" && $inputBand ne "") {
	return ( abs($band - $inputBand) <= 1 );
    }

    # The chromosomes are the same, one or more bands 
    # are not defined
    return 1;
}


sub sameChr {
    # Input: chromosome (in any format)
    # Return: 1 if this object is on the same chromosome

    my ( $self, $inputChr ) = @_;
    my $chr = $self->chr();

    # This should not be called if $chr is not defined
    (defined $chr) ||
	carp $self->symbol(), " does not have a defined chromosome";
  
    # Remove common junk
    $chr =~ s/chr//i; $inputChr =~ s/chr//i;
    $chr =~ s/\_//;  $inputChr =~ s/\_//;
    $chr =~ s/\s+//g; $inputChr =~ s/\s+//g;

    # Change to uppercase in case chromosomes are X or Y
    $chr = uc $chr;
    $inputChr = uc $inputChr;

    if ( $inputChr =~ /^(X|Y|\d+)/ ) {
	$inputChr = $1;
    } else {
	carp "Unexpected format for \$inputChr: '$inputChr' sent to $self->symbol{SYMBOL}";
    }

    if ( $chr =~ /^(X|Y|\d+)/ ) {
	$chr = $1;
    } else {
	carp "Unexpected format for chromosome \$chr: '$chr' in object $self->symbol{SYMBOL}, $self->symbol{CHR}";
    }
    
    return ( $chr eq $inputChr );
}



sub displayObject {
    # For debugging.  Print some of the contents to STDOUT
    my $self = shift;
    
    my ( $idRef );

    print $self->symbol();
    if ( defined $self->id() && $self->id() ne "" ) { print "\tID: ", $self->id(); } 
    if ( defined $self->ensembl() && $self->ensembl() ne "" ) { print "\tEnsembl: ", $self->ensembl(); } 
    $idRef = $self->getAllVegaNames();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tVega: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    $idRef = $self->getAllEntrezNames();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tEntrez: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    $idRef = $self->getAllPreviousSymbols();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tPrevious symbols: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    $idRef = $self->getAllAliases();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tAliases: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    if ( defined $self->chr() && $self->chr() ne "" ) { print "\tChr: ", $self->chr(); } 
}


return 1;




#######################################################
# abreviated pod documentation
#     do    pod2html program.pl > doc.html

#     Add the following line after '<body style="background-color: white">'
#     <div id="TOP"></div>
# the following __END__ token is required if the pod documentation is at the end
__END__



=head1 HugoGene.pm


=head1 SYNOPSIS

Used to hold data from Hugo web site.  These objects are made by HugoGeneMethods::makeHugoGeneObjects

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head1 SUBROUTINES 


=head2 new


=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 id

Set/Get routines for unique variable

=over 4
=item Input parameter(s):


=item Return value(s):

Unique ID number (HGNC ID after removing 'HGNC:')

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 symbol


Set/Get routines for unique variable

=over 4
=item Input parameter(s):


=item Return value(s):

HUGO symbol (gene name)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 name


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Gene description (string with white space)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 type


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Locus type ('pseudogene', 'gene with protein product, inferred', etc.)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 chr


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Chromosome assignment (includes band)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 ensembl


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Ensembl name (mapped data supplied by Ensembl)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 omim


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Online Mendelian Inheritance in Man name (mapped data supplied by NCBI)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 ucsc


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

UCSC ID (mapped data supplied by UCSC)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 start


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Nucleotide start on chr (not supplied by Hugo)

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 stop


Set/Get routines for unique variable
 

=over 4
=item Input parameter(s):


=item Return value(s):

Nucleotide end on chr (not supplied by Hugo) 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 addEntrezName


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 addRefSeqName


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 addVegaName


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 addPreviousSymbol


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 addAlias


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 getAllEntrezNames


Set routines for variables stored in hash
 

=over 4
=item Input parameter(s):


=item Return value(s):

ref to hash 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 getAllRefSeqNames



Describe subroutine here. 

=over 4
=item Input parameter(s):


=item Return value(s):

ref to hash 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 getAllVegaNames



 

=over 4
=item Input parameter(s):


=item Return value(s):

ref to hash 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 getAllPreviousSymbols



 

=over 4
=item Input parameter(s):


=item Return value(s):

ref to hash 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 getAllAliases



 

=over 4
=item Input parameter(s):


=item Return value(s):

ref to hash 

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 isPreviousSymbol



=over 4
=item Input parameter(s):

Hugo name

=item Return value(s):

1 if the input name matches one of the previous symbols of this object

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 isAlias

=over 4
=item Input parameter(s):

Hugo name

=item Return value(s):

1 if the input name matches one of the alias symbols of this object

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 hasVega


=over 4
=item Input parameter(s):

Vega gene ID

=item Return value(s):

1 if this object has the Vega gene ID

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 hasEntrez


=over 4
=item Input parameter(s):

Entrez gene ID
=item Return value(s):

1 if this object has the Entrez gene ID

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 hasRefSeq



=over 4
=item Input parameter(s):

RefSeq gene ID

=item Return value(s):

1 if this object has the RefSeq gene ID

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 overlaps

ASSUMES: we are talking about the same chromosome.  Does not check chromosome.

=over 4
=item Input parameter(s):

start, stop, tolerance (optional)

=item Return value(s):

1 if the start, stop overlap the start, stop of this gene.  undef if this gene does not have start, stop

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 sameMajorChrBand

It does not look at the number after the decimal point. 7p14.20 ~ 7p14 ~ 7p14.2  etc. If one band is 'ter' or 'cen', it will not necessarily return the correct answer since the 'ter' band can have a standard name  

=over 4
=item Input parameter(s):

chromosome (in any format)
=item Return value(s):

1 if this object is on the same major chromosome band

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 sameChr



=over 4
=item Input parameter(s):

chromosome (in any format)
=item Return value(s):

1 if this object is on the same chromosome

=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=head2 displayObject

For debugging.  Print some of the contents to STDOUT 

=over 4
=item Input parameter(s):


=item Return value(s):


=back

=for html <P><A href="#TOP">Top of page</A></FONT><FONT CLASS="textsans9"><P>

=for html <HR>
