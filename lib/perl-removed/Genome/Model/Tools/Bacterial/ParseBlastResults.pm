package Genome::Model::Tools::Bacterial::ParseBlastResults;

use strict;
use warnings;

use Genome;
use FileHandle;
use BPlite;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        blast_query => {
            is  => 'Scalar',
            doc => "blast query fasta file(for length calculations)",
        },
        input => {
            is  => 'Scalar',
            doc => "blast results file",
        },
        num_hits => {
            is  => 'Scalar',
            doc => "number of top hits to report",
        },
        percent => {
            is  => 'Scalar',
            doc => "percent identity",
        },
        fol => {
            is  => 'Scalar',
            doc => "fraction of length",
        },
        output => {
            is  => 'Scalar',
            doc => "output file name",
        },
    ],
    has_optional => [

    ],
);

sub help_brief
{
    "for parsing blast results";
}

sub help_detail
{
    return <<EOS
original usage from old script:

USAGE: parse_blast_results_bit.pl -input <blast results file> -output <output file name> -num_hits <number of tophits to report> -percent <% identity cutoff eg 40> -fol <fraction of query length that must meet % match cutoff> -query <fasta query file (for length calculations)>

NOTE: -fol should be a fractional value, such as .5 for 50% of query length covered, or .25 for 25% or etc...

NOTE: This script only works for queries with standard GenBank names, and only for blast output NOT using the -gi switch


EOS

}

sub help_synopsis
{
    return <<EOS
gmt bacterial parse-blast-results --input ... --output ...
     --num-hits ... --percent .... --fol ... --query-file ...
EOS
}

sub execute
{
    my $self = shift;

    my $query = new FileHandle $self->blast_query;
    my %query_sizes;
    my $current_name;
    while (<$query>)
    {
        chomp;
        my $line = $_;
        next if ( $line =~ /^\s*$/ );
        if ( $line =~ /^>/ )
        {
            ####$line =~ s/^>gi\|\d+\|//; #Get rid of the gi # part of name, since blast (without -gi switch) will not have it either
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            $current_name = $line;
            $current_name =~ s/^>//;
            my @current_name = split( /\s+/, $current_name );
            $current_name = $current_name[0];
        }
        else
        {
            $query_sizes{$current_name} += length($line);
        }
    }

    #Parse blast file
    my $data_hash;
    my $blast           = new FileHandle( $self->input );
    my $multiple_report = new BPlite::Multi( \*$blast );
    while ( my $report = $multiple_report->nextReport )
    {
        my $name = $report->query;
        while ( my $subject = $report->nextSbjct )
        {
            my $subject_name = $subject->name;
            my $best_length  = 0;
            while ( my $hsp = $subject->nextHSP )
            {

#           Collect data in %data_hash ----> NOTE: WE'RE ONLY GRABBING THE LONGEST HSP HERE IF MORE THAN ONE
#           But for Makedonka's purposes it shouldn't matter
                my $current_length = $hsp->length;
                if ( $current_length > $best_length )
                {
                    $best_length = $current_length;
                    $data_hash->{$name}->{$subject_name}->{pvalue} = $hsp->P;
                    $data_hash->{$name}->{$subject_name}->{match}
                        = $hsp->percent;
                    $data_hash->{$name}->{$subject_name}->{sstart}
                        = $hsp->sbjctBegin;
                    $data_hash->{$name}->{$subject_name}->{send}
                        = $hsp->sbjctEnd;
                    $data_hash->{$name}->{$subject_name}->{qstart}
                        = $hsp->queryBegin;
                    $data_hash->{$name}->{$subject_name}->{qend}
                        = $hsp->queryEnd;
                    $data_hash->{$name}->{$subject_name}->{bitscore}
                        = $hsp->bits;
                }
            }
        }
    }
    $blast->close;

    #Prepare output file
    my $output = new FileHandle( ">>" . $self->output );

#Sort data_hash and output the number of top hits requested by user (with min e-value constraint)
    my $i;
    print $output
        "<p_value>  <bitscore>  <%match>  <query start>  <query end>  <subject start>  <subject end>  <subject name>\n";
    print $output
        "(NOTE:displaying only values for longest HSP in each subject)\n";
    foreach my $hit ( keys( %{$data_hash} ) )
    {
        $i = 0;
        foreach my $sbjct (
            sort
            {
                $data_hash->{$hit}->{$a}->{pvalue} <=> $data_hash->{$hit}
                    ->{$b}->{pvalue}
            } keys( %{ $data_hash->{$hit} } )
            )
        {
            if ( $i < $self->num_hits )
            {
                $i++;
                my $pvalue = $data_hash->{$hit}->{$sbjct}->{pvalue};
                my $match  = $data_hash->{$hit}->{$sbjct}->{match};
                my $qstart = $data_hash->{$hit}->{$sbjct}->{qstart};
                my $qend   = $data_hash->{$hit}->{$sbjct}->{qend};
                my $sstart = $data_hash->{$hit}->{$sbjct}->{sstart};
                my $send   = $data_hash->{$hit}->{$sbjct}->{send};
                my $bits   = $data_hash->{$hit}->{$sbjct}->{bitscore};
                my $query_aligned_length = abs( $qend - $qstart );
                my $clean_hit            = $hit;

                if ( $clean_hit =~ /WARNING/ )
                {
                    my $stupidvar = 1;
                }

                my @clean_hit = split( /\(/, $clean_hit );
                $clean_hit = shift(@clean_hit);

                #$clean_hit =~ s/\(\d+\s+letters\;\s+record\s+\d+\)//;
                #$clean_hit =~ s/\(\d+\s+letters\)//;
                $clean_hit =~ s/^\s+//;
                $clean_hit =~ s/\s+$//;
                my @c_h = split( /\s+/, $clean_hit );
                $clean_hit = $c_h[0];

                unless (( defined $query_sizes{$clean_hit} )
                    and ( defined $self->fol ) )
                {
                    my $stupidvar = 1;
                }

                my $min_needed_alignment_length
                    = int( ( $query_sizes{$clean_hit} * $self->fol ) + .5 );
                $sbjct =~ s/^>//;
                if ( $query_aligned_length >= $min_needed_alignment_length )
                {
                    if ( $self->percent <= $match )
                    {
                        if ( $i == 1 )
                        {
                            print $output "========== $hit ==========\n";
                        }
                        print $output
                            "$pvalue\t$bits\t$match\t$qstart\t$qend\t$sstart\t$send\t$sbjct\n";
                    }
                }
            }
        }
    }
    $output->close;

    return 1;
}

1;
