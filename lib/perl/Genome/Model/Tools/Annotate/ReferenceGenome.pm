package Genome::Model::Tools::Annotate::ReferenceGenome;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

my $DEFAULT_OUTPUT_FORMAT = 'gtf';
my $DEFAULT_VERSION = '54_36p_v2';
my $DEFAULT_ANNO_DB = 'NCBI-human.combined-annotation';
my $DEFAULT_REFERENCE_BUILD_ID = 101947881;
my $DEFAULT_SPECIES = 'Homo sapiens';
my $DEFAULT_PICARD_VERSION = '1.36';

class Genome::Model::Tools::Annotate::ReferenceGenome {
    is => ['Command'],
    has => [
        output_file => {
            doc => 'The output file where the annotation file will be dumped.',
        },
    ],
    has_optional => [
        reference_build_id => {
            is => 'Text',
            doc => 'The id of the reference genome build.',
            default_value => $DEFAULT_REFERENCE_BUILD_ID,
        },
        anno_db => {
            doc => 'The name of the annotation database to use. default_value='. $DEFAULT_ANNO_DB,
            default_value => $DEFAULT_ANNO_DB,
        },
        version => {
            doc => 'The version of the annotation database. default_value='. $DEFAULT_VERSION,
            default_value => $DEFAULT_VERSION,
        },
        output_format => {
            doc => 'The file format to output annotation in.',
            valid_values => ['gtf','gff','bed','gff3'],
            default_value => $DEFAULT_OUTPUT_FORMAT,
        },
        species => {
            doc => 'The species we are working with.  Required to get the sequence dictionary.',
            default_value => $DEFAULT_SPECIES,
        },
        picard_version => {
            doc => 'The picard version to create the sequence dictionary if one does not exist.',
            default_value => $DEFAULT_PICARD_VERSION,
        },
    ],
};

sub execute {
    my $self = shift;

    my $output_file = $self->output_file;
    my $dirname = File::Basename::dirname($output_file);
    unless ($dirname) {
        die('Failed to determine output directory from output file: '. $output_file);
    }
    my $tmp_dir = File::Temp::tempdir('Annotate-ReferenceGenome-'.Genome::Sys->username.'-XXXX',DIR => $dirname,CLEANUP => 1);

    my $reference_build = Genome::Model::Build->get($self->reference_build_id);
    unless ($reference_build) {
        die('Failed to get reference build for id: '. $self->reference_build_id);
    }
    my $seq_dict = $reference_build->get_sequence_dictionary('sam',$self->species,$self->picard_version);

    my $tmp_file = Genome::Sys->create_temp_file_path;
    # This is only required to run perl5.10.1 or greater required by Bio-SamTools

    my $cmd = Genome::Model::Tools::BioSamtools::ListChromosomes->execute(
        input_file => $seq_dict,
        output_file => $tmp_file);

    my @chromosomes;
    my $fh = Genome::Sys->open_file_for_reading($tmp_file);
    while (my $line = $fh->getline) {
        chomp($line);
        push @chromosomes, $line;
    }
    $fh->close;
    my %params = (
        anno_db => $self->anno_db,
        version => $self->version,
        chromosomes => \@chromosomes,
        output_format => $self->output_format,
        output_file => $output_file,
        output_directory => $tmp_dir,
    );
    my $module_path = $self->__meta__->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;
    my $workflow = Workflow::Operation->create_from_xml($xml_path);
    Genome::Sys->create_directory($dirname."/annotate_reference_genome_logs");
    $workflow->log_dir($dirname."/annotate_reference_genome_logs");
    my @errors = $workflow->validate;
    unless ($workflow->is_valid) {
        die('Errors encountered while validating workflow '. $xml_path ."\n". join("\n", @errors));
    }
    my $output = Workflow::Simple::run_workflow_lsf($workflow, %params);
    unless (defined $output) {
        @errors = @Workflow::Simple::ERROR;
        for (@errors) {
            print STDERR $_->error ."\n";
        }
        return;
    }
    return 1;
}


1;
