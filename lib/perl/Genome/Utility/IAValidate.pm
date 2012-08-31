package Genome::Utility::IAValidate;

#:eclark 11/17/2009 Code review.

# What is an IA and how do these functions validate it.  The code should be more self-documenting, or have some documentation attached

use strict;
use warnings;
use Text::CSV_XS;
use List::MoreUtils qw/ any /;
use File::Slurp;

sub chromosome
{
    qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
       Y MT NT_113870 NT_113871 NT_113872 NT_113874 NT_113878
       NT_113880 NT_113881 NT_113884 NT_113885 NT_113886 NT_113888
       NT_113889 NT_113890 NT_113898 NT_113899 NT_113901 NT_113902
       NT_113903 NT_113906 NT_113908 NT_113909
       NT_113910 NT_113911 NT_113912 NT_113915 NT_113916
       NT_113917 NT_113923 NT_113924 NT_113925 NT_113926
       NT_113927 NT_113929 NT_113930 NT_113931 NT_113932
       NT_113933 NT_113934 NT_113935 NT_113936 NT_113937 NT_113939
       NT_113943 NT_113944 NT_113946 NT_113949 NT_113951
       NT_113953 NT_113954 NT_113956 NT_113957 NT_113958
       NT_113960 NT_113961 NT_113962 NT_113963 NT_113964
       NT_113965 NT_113966 c5_H2 c6_COX c6_QBL);
}

sub strand { qw(+1 0 -1); }

sub phase { qw(0 1 2); }


sub check_transcript
{
    my $line = shift;
    my $c = Text::CSV_XS->new( { sep_char => "\t" } );
    $c->parse($line);

    my ($tsid,   $gid,    $start,  $stop, $name,
        $source, $status, $strand, $chrom
        )
        = $c->fields();

    # check for dain bramage
    my $flag = 0;
    my @msgs;
    if ( $start == $stop )
    {
        $flag = 1;
        push( @msgs, "start/stop same," );
    }

    unless ( any { $_ eq $chrom; } chromosome() )
    {
        $flag = 1;
        push( @msgs, "chrom out of order, $chrom" );
    }

    unless ( any { $_ eq $strand; } strand() )
    {
        $flag = 1;
        push( @msgs, "strand out of order, $strand" );
    }
    if ( $#msgs > -1 )
    {
        my $msg_string = $line.":".join( ":", @msgs );
        #print STDERR $msg_string,"\n";
        return $msg_string;
    }

    return 1;
}


sub check_sub_structure
{
    my $line = shift;
    my $c = Text::CSV_XS->new( { sep_char => "\t" } );
    $c->parse($line);
    my ( $tid, $sid, $struct, $start, $stop, $ord, $phase, $seq )
        = $c->fields();
    my $flag = 0;
    my @msgs;

    if ( $struct eq 'cds_exon' )
    {

        #carp "phase out of order, $phase\n\t$line"
        #    unless any { $_ eq $phase; } phase();
        unless ( any { $_ eq $phase; } phase() )
        {
            push( @msgs, "phase out of order,$phase" );
            $flag = 1;
        }

        #carp "sequence blank,\n\t$line"
        #    unless ($seq ne '');
        unless ( ( $seq ne '' ) )
        {
            push( @msgs, "sequence blank" );
            $flag = 1;
        }
    }

    if ( $start < 0 )
    {
        push(@msgs, "start position is negative");
    }

    if ( $#msgs > -1 )
    {
        my $msg_string = $line.":".join( ":", @msgs );
        #print STDERR $msg_string;
        return $msg_string;
    }

    return 1;
}


sub check_genes
{
    my $line = shift;
    my $c = Text::CSV_XS->new( { sep_char => "\t" } );
    $c->parse($line);
    my ( $gid, $hugo, $strand ) = $c->fields();
    my @msgs;
    my $flag = 1;
    unless ( any { $_ eq $strand; } strand() )
    {
        push( @msgs, "strand out of order, $strand" );
        $flag = 1;
    }

    if ( $hugo eq '' )
    {
        push( @msgs, "blank hugo name" );
        $flag = 1;
    }

    if ( $#msgs > -1 )
    {
        my $msg_string = $line.":".join( ":", @msgs );
        #print STDERR $msg_string;
        return $msg_string;
    }

    return 1;
}

sub check_protein
{
    my $line = shift;
    my $c = Text::CSV_XS->new( { sep_char => "\t" } );
    $c->parse($line);
    my ($pid, $tid, $name, $seq) = $c->fields();
    my @msgs;
    my $flag = 0;
    my $retval = 1;
    if($seq eq '')
    {
        push(@msgs,"no peptide sequence, $name");
        $flag = 1;
    }
    if( $#msgs > -1 )
    {
        my $msg_string = $line.":".join( ":", @msgs );
        #print STDERR $msg_string,"\n";
        #return $msg_string;
        $retval = $msg_string;
    }
    return $retval;
}

sub check_variant
{
    my $line = shift;

    return 1;
}

sub check_line
{
    my $line = shift;
    my $type = shift;

    if($type =~ /^gene/)
    {
        return check_genes($line);
    }
    elsif($type =~ /^transcript/)
    {
        return check_transcript($line);
    }
    elsif($type =~ /^sub_structure/)
    {
        return check_sub_structure($line);
    }
    elsif($type =~ /^external/)
    {

    }
    elsif($type =~ /^protein/)
    {
        return check_protein($line);
    }
    elsif($type =~ /^variant/)
    {
        return check_variant($line);
    }


    return 1;
}


1;

# $Id$
