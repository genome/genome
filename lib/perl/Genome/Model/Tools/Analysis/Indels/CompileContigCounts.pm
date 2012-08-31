package Genome::Model::Tools::Analysis::Indels::CompileContigCounts;

#This is mostly ported from CompileBowtieResults

use strict;
use warnings;

use IO::File;
use IO::Handle;
use POSIX;    #gonna use SAMTools so will need to check the OS
use Genome;

my %stats = ();

class Genome::Model::Tools::Analysis::Indels::CompileContigCounts {
    is => 'Command',

    has => [
        variant_file    => { is => 'Text', doc => "Indels in annotation format" },
        alignment_files => { is => 'Text', doc => "Alignments in Bam format [comma-separated]" },
        output_file     => { is => 'Text', doc => "Output of indels with alignment counts", is_optional => 1 },
        max_q2_run      => { is => 'Text', doc => "Maximum number of consecutive Q=2 bases before read is discarded", is_optional => 1, default => 30 },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Count valid alignments to indel contigs";
}

sub help_synopsis {
    return <<EOS
This command counts alignments to indel contigs from a bam file
EXAMPLE:	gmt analysis indels compile-contig-counts --variant-file [indels.formatted.tsv] --alignment-files all_seqeunces.bam --output-file [counts.tsv]
EOS
}

sub execute {
    my $self = shift;

    #test that we are on a 64bit system so we can run the GC's installation of samtools
    unless ( POSIX::uname =~ /64/ ) {
        $self->error_message("Must run on a 64 bit machine");
        return;
    }

    ## Get required parameters ##
    my $variant_file    = $self->variant_file;
    my $alignment_files = $self->alignment_files;
    my $output_file     = $self->output_file;
    my $max_q2_run      = $self->max_q2_run;

    $stats{'num_indels'} = $stats{'covered_0x'} = $stats{'covered_1x'} = $stats{'covered_10x'} = $stats{'covered_50x'} = $stats{'covered_100x'} = 0;
    $stats{'num_alignments'} = $stats{'num_alignments_q2_fail'} = $stats{'num_good_alignments'} = 0;

    my $output_fh;
    if ($output_file) {
        $output_fh = IO::File->new( $output_file, "w" );
        unless ($output_fh) {
            $self->error_message("Unable to open $output_file for output");
            return;
        }
    }
    else {
        $output_fh = IO::Handle->new_from_fd( fileno(STDOUT), "w" );
        unless ($output_fh) {
            $self->error_message("Unabke to open a handle to stdout");
            return;
        }
    }

    ## Get the alignment files ##

    my %alignment_totals = ();

    my @alignment_files = split( /\s+,\s+/, $alignment_files );

    foreach my $alignment_file (@alignment_files) {
        print "Parsing $alignment_file\n";

        ## Parse file for alignment counts ##
        my %alignment_counts = parse_alignments( $self, $alignment_file );

        ## Add these counts to total ##

        foreach my $contig ( keys %alignment_counts ) {
            $alignment_totals{$contig} = 0 if ( !$alignment_totals{$contig} );
            $alignment_totals{$contig} += $alignment_counts{$contig};
        }
    }

    ## Parse the variants file ##

    my $input = IO::File->new($variant_file);
    unless ($input) {
        $self->error_message("Unable to open $variant_file");
        return;
    }

    my $lineCounter = 0;

    while (<$input>) {
        chomp;
        my $line = $_;
        $lineCounter++;

        ( my $chrom, my $chr_start, my $chr_stop, my $ref, my $var ) = split( /\t/, $line );
        $chrom     =~ s/[^0-9XYMNT\_]//g;
        $chr_start =~ s/[^0-9]//g if ($chr_start);
        $chr_stop  =~ s/[^0-9]//g if ($chr_stop);

        if ( $chrom && $chr_start && $chr_stop ) {
            $stats{'num_indels'}++;
            my $indel_name = my $indel_type = my $indel_size = my $allele = "";

            if ( $ref eq "0" || $ref eq "-" ) {
                $indel_type = "Ins";
                $indel_size = length($var);
                $allele     = uc($var);

                ## Build indel name ##
                $indel_name = "$chrom\_$chr_start\_$chr_stop\_$indel_type\_$indel_size\_$allele";
            }
            else {
                $indel_type = "Del";
                $indel_size = length($ref);
                $allele     = uc($ref);

                ## Build indel name ##
                $indel_name = "$chrom\_$chr_start\_$chr_stop\_$indel_type\_$indel_size\_$allele";
            }

            my $ref_contig_name = $indel_name . "_ref";
            my $var_contig_name = $indel_name . "_var";

            my $coverage = my $reads1 = my $reads2 = my $var_freq = 0;

            $reads1 = $alignment_totals{$ref_contig_name} if ( $alignment_totals{$ref_contig_name} );
            $reads2 = $alignment_totals{$var_contig_name} if ( $alignment_totals{$var_contig_name} );
            $coverage = $reads1 + $reads2;

            if ($coverage) {
                $var_freq = $reads2 / $coverage * 100;
                $var_freq = sprintf( "%.2f", $var_freq ) . '%';

                $stats{'covered_1x'}++;
                $stats{'covered_10x'}++  if ( $coverage >= 10 );
                $stats{'covered_50x'}++  if ( $coverage >= 50 );
                $stats{'covered_100x'}++ if ( $coverage >= 100 );
            }
            else {
                $stats{'covered_0x'}++;
            }

            print $output_fh join( "\t", $line, $coverage, $reads1, $reads2, $var_freq ) . "\n";
        }
    }

    close($input);

    print STDERR $stats{'num_alignments'} . " alignments parsed\n";
    print STDERR $stats{'num_alignments_q2_fail'} . " had >= $max_q2_run Q2 bases and were discarded\n";
    print STDERR $stats{'num_good_alignments'} . " good alignments were used\n";

    print STDERR $stats{'num_indels'} . " indels in file\n";
    print STDERR $stats{'covered_0x'} . " had no coverage\n";
    print STDERR $stats{'covered_1x'} . " had at least 1x coverage\n";
    print STDERR $stats{'covered_10x'} . " had at least 10x coverage\n";
    print STDERR $stats{'covered_50x'} . " had at least 50x coverage\n";
    print STDERR $stats{'covered_100x'} . " had at least 100x coverage\n";

    return 1;    # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub parse_alignments {
    my ( $self, $bam_file ) = @_;

    my %alignment_counts = ();

    my $max_q2_run = $self->max_q2_run;

    my $q2_string = '#' x $max_q2_run;

    ## Parse the variants file ##

    open BAM, "samtools view -q 1 $bam_file |" or die "Can't fork a pipe for samtools\n";

    my $lineCounter = 0;

    while (<BAM>) {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split( /\t/, $line );
        my $contig_name  = $lineContents[2];
        my $read_quals   = $lineContents[10];

        $stats{'num_alignments'}++;

        if ( $read_quals =~ $q2_string ) {
            $stats{'num_alignments_q2_fail'}++;
        }
        else {
            $alignment_counts{$contig_name}++;
            $stats{'num_good_alignments'}++;
        }
    }

    close(BAM) or die "Error reading BAM file $bam_file\n";

    return (%alignment_counts);
}

1;

