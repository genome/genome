package Genome::Model::Tools::GenePredictor::Rnammer;

use strict;
use warnings;

use Genome;
use Carp 'confess';
use File::Path 'make_path';
use Bio::Tools::GFF;

class Genome::Model::Tools::GenePredictor::Rnammer {
    is => 'Genome::Model::Tools::GenePredictor',
    has_optional => [
        domain => {
            is => 'Text',
            is_input => 1,
            valid_values => ['archaeal', 'bacterial', 'eukaryotic'],
            default => 'eukaryotic',
        },
        version => {
            is => 'Text',
            is_input => 1,
            default => '1.2.1',
            doc => 'Version of rnammer to use',
        },
        molecule_type => {
            is => 'Text',
            doc => 'Specifies molecule types',
            default => 'tsu,lsu,ssu',
        },
        output_format => {
            is => 'Text',
            default => 'gff',
            valid_values => ['fasta', 'gff', 'xml'],
            doc => 'Format of output file',
        },
        temp_working_dir => {
            is => 'Path',
            default => '/tmp/',
            doc => 'Place for temporary files, cleaned up unless keep is set to true',
        },
        debug => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, debugging information is displayed',
        },
        keep => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, temporary files are kept',
        },
        parallel_execution => { 
            is => 'Boolean',
            default => 0,
            doc => 'If set, rnammer is run in parallel',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-M 12000000 -R 'select[type==LINUX64 && mem>12000] rusage[mem=12000]'",
        }
    ],
};

sub help_brief {
    "Write a set of fasta files for an assembly";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip_execution) {
        $self->status_message("Skip execution flag is set, exiting.");
        return 1;
    }

    my $fasta_file = $self->fasta_file;
    $self->status_message("Running rnammer on sequence in $fasta_file");

    unless (-d $self->raw_output_directory) {
        my $mk_rv = make_path($self->raw_output_directory);
        unless ($mk_rv or -d $self->raw_output_directory) {
            confess "Could not make raw ouput directory at " . $self->raw_output_directory;
        }
    }
    $self->status_message("Raw output being placed in " . $self->raw_output_directory);

    # TODO Logic for this output format needs to be added
    if ($self->output_format ne 'gff') {
        $self->error_message("Only GFF output format is currently supported, sorry!");
        confess $self->error_message;
    }

    my $rnammer_path = $ENV{GENOME_SW} . "/rnammer/rnammer-" . $self->version . "/rnammer";
    confess "No rnammer executable found at $rnammer_path!" unless -e $rnammer_path;

    # Create a list of parameters
    my @params;
    push @params, "-T " . $self->temp_working_dir;
    push @params, "-S " . substr($self->domain, 0, 3);
    push @params, "-m " . $self->molecule_type;
    push @params, "-d " if $self->debug;
    push @params, "-multi " if $self->parallel_execution;
    push @params, "-k " if $self->keep;
   
    # Get output file name and create output param for command
    my ($suffix, $param);
    if ($self->output_format eq 'fasta') {
        $suffix .= ".fa";
        $param = "-f ";
    }
    elsif ($self->output_format eq 'gff') {
        $suffix = ".gff";
        $param = "-gff ";
    }
    elsif ($self->output_format eq 'xml') {
        $suffix = ".xml";
        $param = "-xml ";
    }
    my $output_file_fh = File::Temp->new(
        DIR => $self->raw_output_directory,
        TEMPLATE => "rnammer_raw_output_XXXXXX",
        UNLINK => 0,
        CLEANUP => 0,
        SUFFIX => $suffix,
    );
    my $output_file = $output_file_fh->filename;
    $output_file_fh->close;
    chmod(0666, $output_file);
    push @params, $param . $output_file;

    push @params, $fasta_file;

    # Create and execute command
    my $cmd = join(" ", $rnammer_path, @params);
    $self->status_message("Executing rnammer: $cmd");
    my $rna_rv = system($cmd);
    confess "Trouble executing rnammer!" unless defined $rna_rv and $rna_rv == 0;
    $self->status_message("rnammer successfully executed, now parsing output");

    # TODO Add parsing logic for fasta and xml
    if ($self->output_format eq 'gff') {
        # Version 1 is the only version that correctly parses the last column of output...
        my $gff = Bio::Tools::GFF->new(
            -file => $output_file,
            -gff_version => 1,
        );

        my $feature_counter = 0;
        while (my $feature = $gff->next_feature()) {
            $feature_counter++;
            $self->status_message("Parsing " . $feature->seq_id());
            my $gene_name = join(".", $feature->seq_id(), 'rnammer', $feature_counter);
            my ($description) = $feature->get_tag_values('group');

            my $start = $feature->start();
            my $end = $feature->end();
            ($start, $end) = ($end, $start) if $start > $end;

            my $sequence = $self->get_sequence_by_name($feature->seq_id()); 
            confess "Couldn't get sequence " . $feature->seq_id() unless $sequence;
            my $seq_string = $sequence->subseq($start, $end);

            my $rna_gene = Genome::Prediction::RNAGene->create(
                directory => $self->prediction_directory,
                gene_name => $gene_name,
                source => $feature->source_tag(),
                description => $description,
                start => $start,
                end => $end,
                sequence_name => $feature->seq_id(),
                sequence_string => $seq_string,
                strand => $feature->strand(),
                score => $feature->score(),
            );
            confess "Could not create rna gene object!" unless $rna_gene;
        }
    }

    $self->status_message("Parsing done, getting lock and committing!");
    my @locks = $self->lock_files_for_predictions(qw/ Genome::Prediction::RNAGene /);
    UR::Context->commit;
    $self->release_prediction_locks(@locks);

    $self->status_message("Commit done, locks released, rnammer suceessfully completed!");
    return 1;
}

1;
