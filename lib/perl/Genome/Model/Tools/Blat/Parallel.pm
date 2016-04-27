package Genome::Model::Tools::Blat::Parallel;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Blat::Parallel {
    is => ['Genome::Model::Tools::Blat'],
    has => [
            query_files => {
                is_many => 1,
                is_input => 1,
                is => 'ARRAY',
                doc => 'a query file to blat against a genome subject',
            },
            subject_files => {
                is_many => 1,
                is_input => 1,
                is => 'ARRAY',
                doc => 'a list of .fa files to use as subject in blat',
            },
            output_directory => {
                is => 'String',
                doc => 'A network disk directory that all LSF blades can write intermediate files',
                is_input => 1,
            },
            psl_path => {
                is => 'String',
                doc => 'The final blat alignment file concatanated together from each instance',
                is_input => 1,
            },
        ],
    has_optional => [
        blat_params => {
            is => 'string',
            doc => 'a parameter string that will be passed directly to the blat command',
            is_param => 1,
            is_input => 1,
        },
        blat_output_path => {
            is => 'String',
            doc => 'The final output file that contains stderr/stdout from the alginment jobs',
            is_input => 1,
        },
    ],
};

sub help_brief {
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    if (-s $self->psl_path) {
        $self->error_message('A file exists with size for path '. $self->psl_path);
        die;
    }
    unless ($self->blat_output_path) {
        my $blat_output_path = $self->psl_path;
        $blat_output_path =~ s/\.psl/\.out/;
        $self->blat_output_path($blat_output_path);
    }
    if (-s $self->blat_output_path) {
        $self->error_message('A file exists with size for path '. $self->blat_output_path);
        die;
    }
    return $self;
}

sub execute {
    my $self = shift;
    my @subject_files = $self->subject_files;
    my @query_files = $self->query_files;
    my %params = (
                  query_files => \@query_files,
                  subject_files => \@subject_files,
                  psl_path => $self->psl_path,
                  blat_output_path => $self->blat_output_path,
                  output_directory => $self->output_directory,
              );
    if ($self->blat_params) {
        $params{blat_params} = $self->blat_params;
    } else {
        $params{blat_params} = '';
    }
    my $module_path = $self->__meta__->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;

    my $wf = Genome::WorkflowBuilder::DAG->from_xml_filename($xml_path);
    my $output = $wf->execute(inputs => \%params);
    unless (defined $output) {
        $self->fatal_message('Failed to execute workflow.');
    }
    return 1;
}


1;

