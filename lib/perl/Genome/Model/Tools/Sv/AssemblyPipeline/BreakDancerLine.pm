package Genome::Model::Tools::Sv::AssemblyPipeline::BreakDancerLine;

use warnings;
use strict;
use Genome;
use Carp;

class Genome::Model::Tools::Sv::AssemblyPipeline::BreakDancerLine {
    is => 'Genome::Model::Tools::Sv::AssemblyPipeline'
};

# This can accept any of the following formats:


#ID     CHR1    OUTER_START     INNER_START     CHR2    INNER_END       OUTER_END       TYPE    MINSIZE MAXSIZE SOURCE  SCORES
# 1.1     1       1903309 1903309 1       1903656 1903656 DEL     315     315     normal1 90


#ID	CHR1	OUTER_START	INNER_START	CHR2	INNER_END	OUTER_END	TYPE	ORIENTATION	MINSIZE	MAXSIZE	SOURCE	SCORES
# 1.2	1	3252139	3252139	1	3252226	3252226	DEL	+-	43	43	tumor1	109

# 
#CHR1	POS1	CHR2	POS2	ORI	SIZE	TYPE	HET	wASMSCORE	TRIMMED_CONTIG_SIZE	ALIGNED%	NUM_SEG	NUM_FSUB	NUM_FINDEL	BP_FINDEL	MicroHomology	MicroInsertion	PREFIX	ASMPARM
# 2	135855384(135855591)	7	147065745(147065672)	++	208(208)	CTX(CTX)	het	131	441	84%	2	0	0	0	34	-	2.135855591.7.147065672.CTX.208.++	a50.b50

#CHR1	POS1	CHR2	POS2	ORI	SIZE	TYPE	HET	ASMSCORE	ALNSCORE	PREFIX	wASMSCORE	Size_MicroInsertion	ASMPARM
# 1	121053485(121053547)	20	26235211(26235391)	++	328(328)	CTX(CTX)	het	3	28	1.121053547.20.26235391.CTX.328.++	37	0	a50.b100

# Chr1	Pos1	Ori1	Chr2	Pos2	Ori2	Type	Size	Tumor_Het	Tumor_AsmScore	Tumor_AlnScore	Tumor_prefix	Tumor_wASMSCORE	
# 10	100245	-	18	53692	-	CTX	148	het	5	34	10.100172.18.53584.CTX.148.-+	80	-

#CHR1	OUTER_START	CHR2	OUTER_END	TYPE	MINSIZE	ORIENTATION
# 1	6324059	12	7732413	CTX	160	++

# CHR1        POS1        CHR2        POS2        ORT        TYPE        OTHERS
# 8        30265165        17        7108614        ++        CTX        "SJINF001, Germline/Tumor"

sub new {    
    my ($class, $bdLine) = @_; 
    if ( $bdLine !~ /\t/ && $bdLine =~ /\,/ ) { $bdLine =~ s/\,/\t/g; }
    $bdLine =~ s/\"//g;
    bless {
        BD_LINE => $bdLine
    }, $class;
}

sub outputLine {
    my $self = shift;
    return $self->{BD_LINE};
}

sub wAsmScore {
    my $self = shift;
       
    my $line = $self->outputLine();
    my @cols = split /\s+/, $line;
    my $formatType = $self->bdVersion();
    if ( $formatType == 2 || $formatType == 7 ) {
	return $cols[11];
    } elsif ($formatType == 1 || $formatType == 4) {
	return $cols[12];
    }  elsif ($formatType == 3 ) {
	return $cols[8];
    } 
    return undef;
}

sub eventType {
    # This is the assembly-determined event.
    
    my ( @cols, $line, $event, );
    my $self = shift;  
    my $formatType = $self->bdVersion();
  
    $line = $self->outputLine();
    @cols = split /\s+/, $line;

    if ( $formatType == 4 || $formatType == 7 ) {
	$event = $cols[7];
    } elsif ( $formatType == 5 ) {
	$event = $cols[4];
    } elsif ( $formatType == 6 ) {
	$event = $cols[5];
    } else {
	$event = $cols[6];
    }
    if ( $event =~ /(\w+)\(/ ) { $event = $1; }
    ( $event eq "CTX" || $event eq "DEL" || $event eq "INS" || $event eq "INV" || $event eq "ITX" ) || 
	die "Unexpected event type: '$event' from '$line'";

    return $event;
}
    
sub chromosomesAndBreakpoints {
    # Returns ($chr1, $breakPoint1, $chr2, $breakPoint2)
    # $chr1 < $chr2 and when $chr1 eq $chr2, $breakPoint1 < $breakPoint2

    my $self = shift;
    my ($line, @cols, $chr1, $breakPoint1, $chr2, $breakPoint2, $formatType, );
    $line = $self->outputLine();
    @cols = split /\s+/, $line;
    $formatType = $self->bdVersion();
    
    if ( $formatType == 2 || $formatType == 3 ) {
	$chr1 = $cols[0]; $chr2 = $cols[2];
	if ( $cols[1] =~ /(\d+)\(\d+\)/ ) { $breakPoint1 = $1; } else { die "Unexpected format: '$line'"; }
	if ( $cols[3] =~ /(\d+)\(\d+\)/ ) { $breakPoint2 = $1; } else { die "Unexpected format: '$line'"; }

    } elsif ($formatType == 1 ) {
	$chr1 = $cols[0]; $chr2 = $cols[3];
	$breakPoint1 = $cols[2];
	$breakPoint2 = $cols[4];

    } elsif ($formatType == 4 || $formatType == 7) {
	$chr1 = $cols[1]; $chr2 = $cols[4];
	# Use outer boundaries for breakpoints
	$breakPoint1 = $cols[2];
	$breakPoint2 = $cols[6];
    
    } elsif ( $formatType == 5 || $formatType == 6 ) {
	$chr1 = $cols[0]; $chr2 = $cols[2];	
	$breakPoint1 = $cols[1];
	$breakPoint2 = $cols[3];
    }

    if ( ($chr1 =~ /\d+/ && $chr2 =~ /\d+/ && $chr1 > $chr2) ||
	 ($chr1 =~ /[XY]/ && $chr2 =~ /\d+/)
	  || ($chr1 eq $chr2 && $breakPoint1 > $breakPoint2)
	) {
	($chr1, $chr2) = ($chr2, $chr1);
	($breakPoint1, $breakPoint2) = ($breakPoint2, $breakPoint1);
    } 


    return ($chr1, $breakPoint1, $chr2, $breakPoint2);
}

    
sub prefix {
    
    my $self = shift;
    my $prefix;

    my @cols = split /\s+/, $self->{BD_LINE};
    my $formatType = $self->bdVersion();
    
    if ( $formatType == 2 ) {
	$prefix = $cols[10];
    } elsif ( $formatType == 1 ) {
	$prefix = $cols[11];
    } elsif ( $formatType == 3 ) {
	$prefix = $cols[17];
    } elsif ( $formatType == 4 || $formatType == 5 || $formatType == 6) {
	return undef;
    }
    
    ( $prefix =~ /\w+\.\d+\.\w+\.\d+\.\w{3}\.[+-}{+-]/ ) ||
	die "Unexpected format for prefix '$prefix' in '$self->{BD_LINE}'";

    if ( $prefix =~ /(\S+)\.[\+\-]/ ) { $prefix = $1; } 
    return $prefix;

};

sub Id {
    
    my $self = shift;
    my $id;

    my @cols = split /\s+/, $self->{BD_LINE};
    my $formatType = $self->bdVersion();

    ($formatType == 7 || $formatType == 4) || return undef;
    $id = $cols[0];
    if ($id !~ /\w+\.\d+/) { print "Unexpected format for ID: '$id' in ", $self->outputLine(); die; }

    return $id;
}


sub bdVersion {
    # dies if format is different from any of these

    # Returns 7 for this output (this is like '4' except it does not have orientation
    # ID     CHR1    OUTER_START     INNER_START     CHR2    INNER_END       OUTER_END       TYPE    MINSIZE MAXSIZE SOURCE  SCORES
    # 1.1     1       1903309 1903309 1       1903656 1903656 DEL     315     315     normal1 90

    # Returns 6 for this output
    # CHR1        POS1        CHR2        POS2        ORT        TYPE        OTHERS
    # 8        30265165        17        7108614        ++        CTX        "SJINF001, Germline/Tumor"

    # Returns 5 for this output
    # This is the output expected for primer design and pairoscope
    # CHR1	OUTER_START	CHR2	OUTER_END	TYPE	MINSIZE	ORIENTATION
    # 1	6324059	12	7732413	CTX	160	++

    # Returns 4 for this format
    # ID	CHR1	OUTER_START	INNER_START	CHR2	INNER_END	OUTER_END	TYPE	ORIENTATION	MINSIZE	MAXSIZE	SOURCE	SCORES
    # 1.2	1	3252139	3252139	1	3252226	3252226	DEL	+-	43	43	tumor1	109


    # Returns 3 for this format
    #CHR1	POS1	CHR2	POS2	ORI	SIZE	TYPE	HET	wASMSCORE	TRIMMED_CONTIG_SIZE	ALIGNED%	NUM_SEG	NUM_FSUB	NUM_FINDEL	BP_FINDEL	MicroHomology	MicroInsertion	PREFIX	ASMPARM
    # 2	135855384(135855591)	7	147065745(147065672)	++	208(208)	CTX(CTX)	het	131	441	84%	2	0	0	0	34	-	2.135855591.7.147065672.CTX.208.++	a50.b50

    # This is  version of BreakDancer output -- returns 2 for this format
    # CHR1	POS1	CHR2	POS2	ORI	SIZE	TYPE	HET	ASMSCORE	ALNSCORE	PREFIX	wASMSCORE	Size_MicroInsertion	ASMPARM
    # 1	121053485(121053547)	20	26235211(26235391)	++	328(328)	CTX(CTX)	het	3	28	1.121053547.20.26235391.CTX.328.++	37	0	a50.b100


    # This is older version of BreakDancer output  -- returns 1 for this format
    # Chr1	Pos1	Ori1	Chr2	Pos2	Ori2	Type	Size	Tumor_Het	Tumor_AsmScore	Tumor_AlnScore	Tumor_prefix	Tumor_wASMSCORE	
    # 10	100245	-	18	53692	-	CTX	148	het	5	34	10.100172.18.53584.CTX.148.-+	80	-

    
    my $self = shift;
    my ( $line, @cols, );
 
    $line = $self->outputLine();
    @cols = split /\s+/, $line;

    if ( defined $cols[17] && $cols[17] ne "" && $cols[17] =~ /\w+\.\d+\.\w+\.\d+\.\w{3}\.[+-}{+-]/ && $cols[4] =~ /[\+\-]/ ) {
	return 3;
    } elsif ( $cols[4] =~ /[\+\-]/ && $cols[1] =~ /\)/ ) {
	return 2;
    } elsif ( $cols[4] =~ /[\+\-]/ && $cols[1] =~ /^\d+$/ ) {
	return 6;
    } elsif ( $cols[2] =~ /[\+\-]/ && $cols[5] =~ /[\+\-]/ ) {
	return 1;
    } elsif ( $cols[0] =~ /\w+\.\d+/ && defined $cols[8] && $cols[8] =~ /[\+\-]/ && $cols[2] =~ /^\d+$/ && $cols[3] =~ /^\d+$/ ) {
	return 4;
    } elsif ( $cols[0] =~ /\w+\.\d+/ && defined $cols[8] && $cols[8] !~ /[\+\-]/ ) {
	return 7;
    } elsif ( $cols[6] =~ /[\+\-]/ && $cols[1] =~ /^\d+$/ && $cols[3] =~ /^\d+$/ ) {
	return 5;
    } else {
	die "Unexpected format: '$line'";
    }

}   


sub sameRegion {
    # Gets ref to BreakDancerLine and optional $buffer size (default is 0)
    # Returns 1 if the regions defined by this line and the input line are the
    # same within $buffer

    my $self = shift;
    my ($bdRef, $buffer) = @_;
    if ( !defined $buffer ) { $buffer = 0; }

    # Need to be same event
    ( $self->eventType() eq $bdRef->eventType() ) || return 0;

    my ($chr1, $breakPoint1, $chr2, $breakPoint2) = $bdRef->chromosomesAndBreakpoints();
    my ($thisChr1, $thisBreakPoint1, $thisChr2, $thisBreakPoint2) = $self->chromosomesAndBreakpoints();
    
    return ( $chr1 eq $thisChr1 && $chr2 eq $thisChr2 &&
	     abs($breakPoint1 - $thisBreakPoint1) <= $buffer &&
	     abs($breakPoint2 - $thisBreakPoint2) <= $buffer
	);
}

sub regionsOverlap {
    # Gets ref to BreakDancerLine and optional $buffer size (default is 0)
    # Returns 1 if the regions defined by this line and the input line overlap within $buffer

    my $self = shift;
    my ($bdRef, $buffer) = @_;
    if ( !defined $buffer ) { $buffer = 0; }

    # Need to be same event
    ( $self->eventType() eq $bdRef->eventType() ) || return 0;

    my ($chr1, $breakPoint1, $chr2, $breakPoint2) = $bdRef->chromosomesAndBreakpoints();
    my ($thisChr1, $thisBreakPoint1, $thisChr2, $thisBreakPoint2) = $self->chromosomesAndBreakpoints();
    
    # Same chromosome(s)
    if ( $chr1 ne $thisChr1 || $chr2 ne $thisChr2 ) { return 0; }

    # Make sure first breakpoint < second breakpoint for both before checking for overlap
    if ( $breakPoint1 > $breakPoint2 ) { ($breakPoint1, $breakPoint2) = ($breakPoint2, $breakPoint1); }
    if ( $thisBreakPoint1 > $thisBreakPoint2 ) { ($thisBreakPoint1,$thisBreakPoint2) = ($thisBreakPoint2,$thisBreakPoint1); }
    
    # return ( $stop + $tolerance > $thisStart && $start - $tolerance < $thisStop );

    return ( $breakPoint2 + $buffer > $thisBreakPoint1 && $breakPoint1 - $buffer < $thisBreakPoint2);


}



return 1;
