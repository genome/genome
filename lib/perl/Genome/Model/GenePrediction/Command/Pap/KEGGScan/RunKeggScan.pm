package Genome::Model::GenePrediction::Command::Pap::KEGGScan::RunKeggScan;

# bdericks: This module was originally written by Todd Wylie and named
# KEGGscan_KO.v50.pl. The module can be found in version control at 
# svn+ssh://svn/srv/svn/gscpan/hgmi_annotation/proteinAnnot. I've made
# some small changes so this script can be included as a tool in the PAP
# namespace and included in all further snapshots.

use strict;
use warnings;

use PAP;

use Carp qw(confess);
use Cwd;
use File::Copy;
use File::Basename;
use File::Path qw(make_path);

use Genome::Model::GenePrediction::Command::Pap::KEGGScan::KeggScanCommon
        qw( create_directory
            revise_subject_for_ECs
            build_EC_index
            create_report_line
            generate_blast_reports
        );

# FIXME Really need to remove this dependency, but porting this script
# into the PAP namespace may not be straightforward... Used to generate
# the top and full reports after blast is run
use lib "/gsc/scripts/gsc/compbio/lib";
use BPdeluxe;

class Genome::Model::GenePrediction::Command::Pap::KEGGScan::RunKeggScan {
    is => 'Command::V1',
    has => [
        species_name => {
            is => 'Text',
            doc => 'Name of input species',
        },
        query_fasta_path => {
            is => 'Path',
            doc => 'Path to query fasta file',
        },
        subject_fasta_path => {
            is => 'Path',
            doc => 'Path to subject fasta file',
        },
        output_directory => {
            is => 'Path',
            doc => 'All output is placed in this directory',
        },
    ],
    has_optional => [
        query_sequence_type => {
            is => 'Text',
            default => 'CONTIG',
            doc => 'Type of query sequence',
        },
        kegg_release_version => {
            is => 'Text',
            default => 'RELEASE-52',
            doc => 'Version of Kegg to use',
        },
        blast_lsf_queue => {
            is => 'Text',
            default => $ENV{GENOME_LSF_QUEUE_SHORT},
            doc => 'Queue in which LSF jobs should be scheduled',
        },
        blast_lsf_resource => {
            is => 'Text',
            default => 'select[mem>1024] rusage[mem=1024]',
            doc => 'Resources needed to run LSF jobs',
        },
        blast_lsf_max_memory => {
            is => 'Number',
            default => '1024000',
            doc => 'Maximum memory limit of LSF jobs',
        },
        blast_lsf_job_limit => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of concurrent LSF jobs allowed',
        },
        fasta_chunk_size => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of sequence allowed in a fasta chunk',
        },
        db_format => {
            is => 'Text',
            calculate_from => 'kegg_release_version',
            calculate => q|
                $kegg_release_version =~ /\D*(\d+)\D*/;
                my $version_num = $1;
                return "new" if $version_num >= 41;
                return "old"; |,
            doc => 'The format of KO entries in genes db changed in version 41, ' .
                   'older versions require special handling',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->status_message("Creating kegg output directory at " . $self->output_directory);

    if (-d $self->output_directory) {
        $self->warning_message("Existing output directory found at " . $self->output_directory . ", removing!");
        my $rv = system("rm -rf " . $self->output_directory);
        confess "Could not remove " . $self->output_directory unless defined $rv and $rv == 0;
    }

    make_path($self->output_directory);
    chmod(0775, $self->output_directory);

    $self->status_message("Done creating output directory, now creating EC indices.");

    # We will be running BLAST (blastp) on the set of all sequences in the
    # query (species) FASTA file against ONLY the sequences in the KEGG
    # genes FASTA file that have associated EC numbers. Place the EC-only
    # FASTA file in the session output directory. We will be splitting up
    # the BLAST jobs across multiple CPUs via the blade center and
    # BladeBlastBatcher.pm.
    my $revised_subject = $self->revise_subject_for_ECs;
    my ($EC_index,$KO_index) = $self->build_EC_index($revised_subject);

    $self->status_message("Using BladeBlastBatcher to schedule blast jobs.");

    my $blast_batcher = Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher->create(
	    query_fasta_path => $self->query_fasta_path,
		subject_fasta_path => $revised_subject, 
        lsf_queue => $self->blast_lsf_queue,   
		lsf_job_limit => $self->blast_lsf_job_limit, 
		output_directory => $self->output_directory, 
		blast_name => "blastp",
		blast_params => "hitdist=40 wordmask=seg postsw B=10000 topcomboN=1",
        lsf_resource => $self->blast_lsf_resource,
        lsf_max_memory => $self->blast_lsf_max_memory,
        fasta_chunk_size => $self->fasta_chunk_size,
	);
    my $batcher_rv = $blast_batcher->execute;
    confess "Trouble executing blast batcher!" unless defined $batcher_rv and $batcher_rv == 1;
    my @reports = @{$blast_batcher->reports};

    $self->status_message("BladeBlastBatcher completed, aggregating reports.");
    $self->status_message("There are " . scalar @reports . " that need to be aggregated!");

    my $master_report_path = $self->output_directory . "/MASTER_BLAST.report";
    my $aggregate = IO::File->new($master_report_path, "a");
    REPORT: for my $report (@reports) {
        unless (-e $report) {
            $self->warning_message("Report at $report doesn't exist!");
            next REPORT;
        }

        unless (-s $report) {
            $self->warning_message("Report at $report has no size!");
            unlink $report;
            next REPORT;
        }

        my $fh = new IO::File->new($report, "r");
        while (my $line = $fh->getline) { 
            $aggregate->print($line);
        }
        $fh->close;
        unlink($report);   # Saves a little space.
    }
    $aggregate->close;

    unless (-s $master_report_path) {
        confess "Master blast report at $master_report_path has no size!";
    }

    $self->status_message("All reports aggregated, now parsing.");

    # Walk through all of the blast report files and gather the needed info. We
    # will be using BPdeluxe.pm to parse the BLAST reports. Our goal is to get
    # one best (top) match per query with regard to E-value & bit score. If
    # there are multiple/identical E-value matches under a single query, then
    # best bitscore will determine the top hit. If both E-value and bitscore
    # are identical we will simply use the first match in the list by default.
    my $full_report_path    = $self->output_directory . "/" . "REPORT-full.ks";
    my $top_hit_report_path = $self->output_directory . "/" . "REPORT-top.ks";
    $self->generate_blast_reports($master_report_path, $top_hit_report_path, $full_report_path, $EC_index, $KO_index);
    $self->status_message("Report parsing completed, starting clean up.");

    # Move files to subdirectories to reduce clutter
    my $ancillary_dir = $self->create_directory($self->output_directory . "/ancillary/");
    my $log_dir = $self->create_directory($self->output_directory . "/logs/");

    opendir(DIR, $self->output_directory) or confess "Couldn't open " . $self->output_directory;
    while (my $file = readdir(DIR)) {
        my $destination;
        if ($file =~ /^CHUNK/) {
            unlink ($self->output_directory . "/$file");
        }
        elsif (($file =~ /OU$/) || ($file =~ /EC_only/)) {
            $destination = $ancillary_dir
        }
        elsif (($file =~ /.out$/) || ($file =~ /.err$/) || ($file =~ /.log$/)) {
            $destination = $log_dir;
        }
        next unless defined $destination;

        my $full_path = $self->output_directory . "/$file";
        my $mv_rv = system("mv $full_path $destination");
        $self->warning_message("Could not move $full_path to $destination!") unless $mv_rv == 0;
    }

    $self->status_message("Keggscan finished, reports can be found at:\n$full_report_path\n$top_hit_report_path");
    return 1;
}

sub bpdelux_module_name { q(BPdeluxe::Multi) }

# *|**********************************************************************|*
# *| POD: Section below is reserved for documentation.                    |*
# *| Last update: Thu Oct 20 10:08:34 CDT 2005                            |*
# *|**********************************************************************|*

=head1 NAME

RunKeggScan.pm

=head1 USAGE

 Your query sequence file should be of nucleotide FASTA format. The subject database should be a FASTA of protein sequences as retrieved from the KEGG FTP site (ftp://ftp.genome.ad.jp/pub/kegg/). As of this writing, the genes file can be found by downloading and unpacking sequences.tar.gz from KEGG; the file needed is simply called "genes".

 NOTE: The 'KeggRelease' line entry must ONLY contain a single numeric portion.  This code uses this number to determine the expected formatting of the genes db fasta file with regards to the KO numbers.

=head1 DESCRIPTION

KEGGscan was written to compare nucleotide sequence (full length cDNAs, genes, clustered ESTs) to the KEGG gene database of /currated protein sequences. Sequences in the KEGG database have known, annotated Enzyme Commission (EC) ids associated with them. By aligning query sequence against annotated sequence we can assign putative function by EC number association. Obviously, the better the resolution of the query sequence, the more believable our results will be.  The BLAST used to make our alignments uses the topcomboN=1 argument to return only the "best" result. A routine takes the KEGG genes FASTA file and revises it, creating a version with only entries that have associated Enzyme Commission (EC) numbers.

The following logic filters the match list. Top hits are determined on a one-to-one basis (1 group hit per query). The rules for hits are:

=over 5

=item 1

Lowest p_value determines top match.                                   

=item 2

If other HSPs tie the lowest p_value, highest bit score prevails.      

=item 3

If p_value, bit scores are identical, then choosing any 1 of the hits is equally representative.              

=back

=head1 AUTHOR

 Todd Wylie
 CompBio
 Genome Sequencing Center
 Washington University 
 School of Medicine
 4444 Forest Park Boulevard
 St. Louis, MO 63108

CONTACT: twylie@.wustl.edu

 Modifications made by Brian Derickson
 Analysis Pipeline Group
 The Genome Center at Washington University in St. Louis

CONTACT: bdericks@genome.wustl.edu

=head1 LIMITATION/BUGS

Requires:

=over 5

=item 1

BladeBlastBatcher.pm

=item 2

BlastTopHitLogic.pm  

=back

NOTE: Looks like KEGG genes files always have subject_name in 1st col & EC in last. If this ever changes, we'll have to update this routine....

=head1 COPYRIGHT

Copyright (C) 2005 by Todd Wylie and Washington University School of Medicine 
Genome Sequencing Center.

=head1 NOTE:

This software was written using the latest version of GNU Emacs, the extensible, 
real-time text editor. Please see http://www.gnu.org/software/emacs/ for more 
information and download sources.

=cut
