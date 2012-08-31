package Genome::Model::Tools::SmrtAnalysis::BaseModSummary;

use strict;
use warnings;
use Statistics::Descriptive;
use Genome;

# T. Wylie  <twylie@genome.wustl.edu>
# Mon May 21 14:14:11 CDT 2012

class Genome::Model::Tools::SmrtAnalysis::BaseModSummary {
    is        => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
                  ref_fasta => { doc => 'The reference FASTA used to generate the base mod GFF file.' },
                  gff       => { doc => 'The PacBio base mod GFF file' },
                  mode      => {
                                doc           => 'Run mode for reporting, default is REPORT. [REPORT|PARSE]',
                                default_value => 'REPORT',
                               },
                  minscore  => {
                                doc           => 'Minimum base mod call score required for membership into report, default is 0.',
                                default_value => '0',
                               },
                  prepos    => {
                                doc           => 'Base position to start the motif string before the call position, default is 3.',
                                default_value => '3',
                               },
                  poslen    => {
                                doc           => 'Length of motif string starting from PREPOS value, default is 4.',
                                default_value => '4',
                               },
                 },
};

sub help_brief {
    'Parse and summarize PacBio base modification .gff files.'
}

sub help_detail {
    return <<EOS
Parse and summarize PacBio base modification .gff files. Requires the reference FASTA on which the .gff base mod file was built. Can be run in either REPORT or PARSE modes.
EOS
}

sub execute {
    my $self = shift;

    # INPUT
    my $mode     = $self->mode();
    my $ref      = $self->ref_fasta();
    my $file     = $self->gff();
    my $minscore = $self->minscore();
    my $precall  = $self->prepos();
    my $postcall = $self->poslen();

    # Index the reference FASTA file.
    my $ref_sequence;
    open (REF, "$ref");
    while (<REF>) {
        chomp;
        if ($_ =~ /\>/) { next }
        $ref_sequence .= $_;
    }
    close (REF);

    my (%report, %motif_metrics, %call_metrics);

    # Reading/parsing the base mod GFF file.
    open (IN, "$file");
    LINE:
    while (<IN>) {
        chomp;
        next LINE if ($_ =~ /^\#/);
        my @f             = split( "\t", $_ );
        my $position      = $f[3];
        my $score         = $f[5];
        my $strand        = $f[6];
        my $glob          = $f[8];
        my $refcall       = substr( $ref_sequence, $position - 1, 1 );
        my $refcall_motif;
        my ($IPDRatio, $context, $coverage) = split( /\;/, $glob);
        $coverage =~ s/coverage\=//g;

        if ($strand eq '+') {
            $refcall_motif = substr( $ref_sequence, $position - ($precall + 1), $postcall );
        }
        elsif ($strand eq '-') {
            $refcall_motif = substr( $ref_sequence, $position - ($postcall - $precall), $postcall );
            $refcall_motif = reverse_complement( $refcall_motif );
            $refcall       = reverse_complement( $refcall );
        }

        if ($score >= $minscore) {
            # Parse mode.
            print join(
                       "\t",
                       $position,
                       $strand,
                       $score,
                       $coverage,
                       $refcall,
                       $refcall_motif,
                      ) . "\n" if ($mode =~ /PARSE/i);

            # Report mode.
            $report{'motif'}->{$refcall_motif}++;
            $report{'motif_total'}++;
            $report{'call'}->{$refcall}++;
            $report{'call_total'}++;
            push (@{ $motif_metrics{$refcall_motif}->{'score'}    }, $score);
            push (@{ $motif_metrics{$refcall_motif}->{'coverage'} }, $coverage);
            push (@{ $motif_metrics{$refcall_motif}->{'positions'} }, $position);
            $motif_metrics{$refcall_motif}->{'plus_strand'}++  if ($strand eq '+');
            $motif_metrics{$refcall_motif}->{'minus_strand'}++ if ($strand eq '-');
            push (@{ $call_metrics{$refcall}->{'score'}    }, $score);
            push (@{ $call_metrics{$refcall}->{'coverage'} }, $coverage);
            $call_metrics{$refcall}->{'plus_strand'}++  if ($strand eq '+');
            $call_metrics{$refcall}->{'minus_strand'}++ if ($strand eq '-');

        }
    }
    close (IN);

    if ($mode =~ /REPORT/i) { print_report( \%report, \%motif_metrics, \%call_metrics, $minscore ) }

    return 1;  # end of execute
}

sub round { return sprintf( "%.2f", shift ) };

sub order_array {
    my $aref = shift;
    my %unique;
    my @ordered_array;
    foreach my $pos (@{ $aref }) { $unique{$pos} = $pos }
    foreach my $hpos (sort {$unique{$a} <=> $unique{$b}} keys %unique) {
        push( @ordered_array, $hpos );
    }
    return \@ordered_array;
}

sub reverse_complement {
    my $string = shift;
    my @letters = split( '', reverse( $string ));
    my $reverse_complement;
    foreach my $letter (@letters) {
        if ($letter =~ /A/i) {
            $reverse_complement .= 'T';
        }
        elsif ($letter =~ /C/i) {
            $reverse_complement .= 'G';
        }
        elsif ($letter =~ /G/i) {
            $reverse_complement .= 'C';
        }
        elsif ($letter =~ /T/i) {
            $reverse_complement .= 'A';
        }
    }
    return $reverse_complement;
}

sub print_report {
    my ($report, $motif_metrics, $call_metrics, $minscore) = @_;

    my %report        = %{ $report };
    my %motif_metrics = %{ $motif_metrics };
    my %call_metrics  = %{ $call_metrics };

    # CALLS
    my $total_calls = keys( %{ $report{'call'} } );
    print "**BASE MOD CALLS**\n\n";
    print join(
               "\t",
               'MINSCORE',
               'CALL',
               'CALL COUNT',
               'TOTAL COUNT',
               '% TOTAL',
               'MINUS STRAND',
               '% MINUS STRAND',
               'PLUS STRAND',
               '% PLUS STRAND',
               'AVE COV',
               'STDEV COV',
               'AVE SCORE',
               'STDEV SCORE',
              ) . "\n";
    foreach my $call (sort {$report{'call'}->{$b} <=> $report{'call'}->{$a}} keys %{ $report{'call'} }) {
        my $percentage   = ($report{'call'}->{$call} / $report{'call_total'}) * 100;

        my $myCallStats   = Statistics::Descriptive::Full->new();
        $myCallStats->add_data( @{ $call_metrics{$call}->{'coverage'} } );
        my $coverage_ave   = round( $myCallStats->mean() );
        my $coverage_stdev = round( $myCallStats->standard_deviation() );

        $myCallStats    = Statistics::Descriptive::Full->new();
        $myCallStats->add_data( @{ $call_metrics{$call}->{'coverage'} } );
        my $score_ave   = round( $myCallStats->mean() );
        my $score_stdev = round( $myCallStats->standard_deviation() );

        if (!$call_metrics{$call}->{'minus_strand'}) { $call_metrics{$call}->{'minus_strand'} = 0 }
        if (!$call_metrics{$call}->{'plus_strand'})  { $call_metrics{$call}->{'plus_strand'}  = 0 }
        my $total_strand     = $call_metrics{$call}->{'minus_strand'} + $call_metrics{$call}->{'plus_strand'};
        my $plus_percentage  = round( ($call_metrics{$call}->{'plus_strand'}  / $total_strand) * 100 ) . '%';
        my $minus_percentage = round( ($call_metrics{$call}->{'minus_strand'} / $total_strand) * 100 ) . '%';

        print join(
                   "\t",
                   '@' . $minscore,
                   $call,
                   $report{'call'}->{$call},
                   $report{'call_total'},
                   round( $percentage ) . '%',
                   $call_metrics{$call}->{'minus_strand'},
                   $minus_percentage,
                   $call_metrics{$call}->{'plus_strand'},
                   $plus_percentage,
                   $coverage_ave,
                   $coverage_stdev,
                   $score_ave,
                   $score_stdev,
                  ) . "\n";
    }
    print "\n\n";

    # MOTIFS
    my $total_motifs = keys( %{ $report{'motif'} } );
    print "**BASE MOD MOTIFS**\n";
    print "\n";
    print join(
               "\t",
               'MINSCORE',
               'MOTIF',
               'MOTIF COUNT',
               'TOTAL COUNT',
               '% TOTAL',
               'MINUS STRAND',
               '% MINUS STRAND',
               'PLUS STRAND',
               '% PLUS STRAND',
               'AVE COV',
               'STDEV COV',
               'AVE SCORE',
               'STDEV SCORE',
               'POSITIONS',
              ) . "\n";
    foreach my $motif (sort {$report{'motif'}->{$b} <=> $report{'motif'}->{$a}} keys %{ $report{'motif'} }) {
        my $percentage = ($report{'motif'}->{$motif} / $report{'motif_total'}) * 100;

        my $myMotifStats   = Statistics::Descriptive::Full->new();
        $myMotifStats->add_data( @{ $motif_metrics{$motif}->{'coverage'} } );
        my $coverage_ave   = round( $myMotifStats->mean() );
        my $coverage_stdev = round( $myMotifStats->standard_deviation() );

        $myMotifStats      = Statistics::Descriptive::Full->new();
        $myMotifStats->add_data( @{ $motif_metrics{$motif}->{'score'} } );
        my $score_ave   = round( $myMotifStats->mean() );
        my $score_stdev = round( $myMotifStats->standard_deviation() );

        if (!$motif_metrics{$motif}->{'minus_strand'}) { $motif_metrics{$motif}->{'minus_strand'} = 0 }
        if (!$motif_metrics{$motif}->{'plus_strand'})  { $motif_metrics{$motif}->{'plus_strand'}  = 0 }
        my $total_strand     = $motif_metrics{$motif}->{'minus_strand'} + $motif_metrics{$motif}->{'plus_strand'};
        my $plus_percentage  = round( ($motif_metrics{$motif}->{'plus_strand'}  / $total_strand) * 100 ) . '%';
        my $minus_percentage = round( ($motif_metrics{$motif}->{'minus_strand'} / $total_strand) * 100 ) . '%';

        my @positions = @{ order_array( $motif_metrics{$motif}->{'positions'} ) };

        print join(
                   "\t",
                   '@' . $minscore,
                   $motif,
                   $report{'motif'}->{$motif},
                   $report{'motif_total'},
                   round( $percentage ) . '%',
                   $motif_metrics{$motif}->{'minus_strand'},
                   $minus_percentage,
                   $motif_metrics{$motif}->{'plus_strand'},
                   $plus_percentage,
                   $coverage_ave,
                   $coverage_stdev,
                   $score_ave,
                   $score_stdev,
                   join( ";", @positions ),
                  ) . "\n";
    }
    print "\n\n";

}

1;  # end of package


__END__

**EXAMPLE FILE**

##gff-version 3
##source kineticModificationDetector 0.1
##source-commandline /gscmnt/pacbio/production/smrtanalysis-test//analysis/bin/ipdSummary.py --control /gscmnt/pacbio/production/smrtanalysis-test/common/jobs/016/016454/data/aligned_reads.cmp.h5 --summary_h5 /gscmnt/pacbio/production/smrtanalysis-test/common/jobs/016/016456/data/temp_kinetics.h5 --gff /tmp/tmpi6PcU7.gff --csv /tmp/tmpQWUZ1Q.csv --reference /gscmnt/pacbio/production/smrtanalysis-test/common/references/lambda /gscmnt/pacbio/production/smrtanalysis-test/common/jobs/016/016456/data/aligned_reads.cmp.h5
##sequence-header ref000001 lambda_NEB3011
    ref000001kinModCallmodified_base116801168020+.IPDRatio=1.24;context=GCTGAGCAGCAGACTCAACAGGACAAAAATGCGCAGCAGCA;coverage=353.00
    ref000001kinModCallmodified_base119541195420-.IPDRatio=1.20;context=AAGCGTCAGCAGGGCAGCATGAGCACTGTCTTCCTGACGAT;coverage=330.00
...
