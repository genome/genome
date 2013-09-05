package Genome::Model::Tools::Predictor::Keggscan;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use File::Path qw(make_path);

#FIXME This file has several packages declared inside
#It was copied from /gsc/scripts/gsc/compbio/lib to
#eliminate the untracked external dependency and
#prefixes to the package name were added.
use Genome::Model::Tools::Predictor::BPdeluxe; 

use Genome::Model::GenePrediction::Command::Pap::KEGGScan::KeggScanCommon
        qw( create_directory
            revise_subject_for_ECs
            build_EC_index
            create_report_line
            generate_blast_reports
        );

class Genome::Model::Tools::Predictor::Keggscan {
    is => 'Genome::Model::Tools::Predictor::Base',
    doc => 'execute keggscan',

    #parameters originally from RunKeggscan but called here with constant values
    has_constant => [
        species_name => {
            is => 'Text',
            doc => 'Name of input species',
            value => 'default',
        },
        blast_lsf_queue => {
            is => 'Text',
            value => 'long',
            doc => 'Queue in which LSF jobs should be scheduled',
        },
        blast_lsf_resource => {
            is => 'Text',
            value => 'select[mem>4096] rusage[mem=4096]',
            doc => 'Resources needed to run LSF jobs',
        },
        blast_lsf_max_memory => {
            is => 'Number',
            value => '4096000',
            doc => 'Maximum memory limit of LSF jobs',
        },
        blast_lsf_job_limit => {
            is => 'Number',
            value => 50,
            doc => 'Maximum number of concurrent LSF jobs allowed',
        },
        fasta_chunk_size => {
            is => 'Number',
            default => 10,
            doc => 'Maximum number of sequence allowed in a fasta chunk',
        },
        kegg_release_version => {
            is => 'Text',
            value => 'RELEASE-52',
            doc => 'Version of Kegg to use',
        },
        query_sequence_type => {
            is => 'Text',
            value => 'CONTIG',
            doc => 'Type of query sequence',
        },
    ],
    has_calculated => [
        subject_fasta_path => {
            is => 'Path',
            doc => 'Path to subject fasta file',
            calculate_from => ['version'],
            calculate => q{
                # FIXME Figure out what this is... it needs to be removed.
                return "/gscmnt/temp212/info/annotation/KEGG/Version_$version/genes.v$version.faa";
            },
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

sub requires_chunking {
    return 0;
}

# tmooney: The bulk of this was previously a separate module,
# `gmt predictor keggscan run-keggscan`, but is now incorporated here.
# The following message was in that module.
# bdericks: This module was originally written by Todd Wylie and named
# KEGGscan_KO.v50.pl. The module can be found in version control at 
# svn+ssh://svn/srv/svn/gscpan/hgmi_annotation/proteinAnnot. I've made
# some small changes so this script can be included as a tool in the PAP
# namespace and included in all further snapshots.

sub run_predictor {
    my $self = shift;

    $self->status_message("Starting run");
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

    my $blast_batcher = Genome::Model::Tools::Predictor::Keggscan::BladeBlastBatcher->create(
        query_fasta_path => $self->input_fasta_file,
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

sub parse_output {
    my $self = shift;

    my $top_output_fh = IO::File->new($self->raw_output_path);
    unless ($top_output_fh) {
        die "Could not open file handle for raw output file " . $self->raw_output_path;
    }
    
    my @features;
    LINE: while (my $line = $top_output_fh->getline) {
        chomp $line;

        my @fields = split /\t/, $line;
        my (
            $ec_number,
            $gene_name,
            $hit_name,
            $e_value,
            $description,
            $orthology
        ) = @fields[1,2,3,4,6,8];

        # The values in the third column should be in this form:  gene_name (N letters; record M).
        ($gene_name) = split /\s+/, $gene_name; 

        # Some descriptions have EC numbers embedded in them, which should be removed
        $description =~ s/\(EC .+\)//;
        
        my $feature = new Bio::SeqFeature::Generic(-display_name => $gene_name);
        $feature->add_tag_value('kegg_evalue', $e_value);
        $feature->add_tag_value('kegg_description', $description);
        
        my $gene_dblink = Bio::Annotation::DBLink->new(
            -database   => 'KEGG',
            -primary_id => $hit_name,
        );
        $feature->annotation->add_Annotation('dblink', $gene_dblink);

        # Sometimes there is no orthology identifier (value is literally 'none').
        # It is not unforseeable that it might also be blank/undefined.  
        if (defined($orthology) && ($orthology ne 'none')) {
            my $orthology_dblink = Bio::Annotation::DBLink->new(
                -database   => 'KEGG',
                -primary_id => $orthology,
            );
            $feature->annotation->add_Annotation('dblink', $orthology_dblink);
        }
        
        push @features, $feature;
    }
    
    $self->unfiltered_bio_seq_features(\@features);
    return 1;
}

# This was taken from the kegg2ace, which uses the raw output of keggscan. This should
# be refactored to operate directly on the Bio::SeqFeature objects stored in the
# filtered_bio_seq_features property. For now, the filter logic is duplicated here.
# Thankfully this logic is really simple. :)
sub create_ace_file {
    my $self = shift;
    my $fh = Genome::Sys->open_file_for_reading($self->raw_output_path);
    my $ace_fh = Genome::Sys->open_file_for_writing($self->ace_file_path);

    while (my $line = $fh->getline) {
        chomp $line;
        my @fields = split(/\t/, $line);
        
        my $gene = $fields[2];
        $gene =~ s/\(.*//;
        my $code = $fields[3];
        my $e_value = $fields[4];
        my $desc = $fields[6];
        $desc =~ s/\(EC .+\)//;
        my $ko = $fields[8];

        next if $e_value >= 0.01;

        if ($ko eq 'none') {
            $ko = "";
        }

        $ace_fh->print("Sequence \"$gene\"\nKEGG   \"$code $e_value $desc $ko\"\n\n");
    }

    return 1;
}

sub filter_results {
    my $self = shift;
    my @features;
    for my $feature ($self->unfiltered_bio_seq_features) {
        my ($evalue) = $feature->get_tag_values('kegg_evalue');
        next if $evalue >= 0.01;
        push @features, $feature;
    }
    $self->bio_seq_features(\@features);
    return 1;
}

sub tool_path_for_version {
    my ($self, $version) = @_;
}

sub raw_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'REPORT-top.ks');
}

sub debug_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'keggscan.debug');
}

sub dump_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'keggscan.output');
}

sub ace_file_path {
    my $self = shift;
    return join('/', $self->output_directory, 'keggscan.ace');
}

sub bpdelux_module_name { q(Genome::Model::Tools::Predictor::BPdeluxe::Multi) }

1;

