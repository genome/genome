package Genome::Model::Tools::Blat::Subjects;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Blat::Subjects {
    is => ['Genome::Model::Tools::Blat'],
    has => [
            query_file => {
                           is_input => 1,
                           is => 'String',
                           doc => 'a query file to blat against a genome subject',
                      },
            subject_files => {
                              is_many => 1,
                              is_input => 1,
                              is => 'String',
                              doc => 'a list of .fa files to use as subject in blat',
                          },
            psl_path => {
                         is => 'String',
                         is_input => 1,
                     },
        ],
    has_optional => [
                     blat_params => {
                                     is => 'string',
                                 },
                     blat_output_path => {
                                          is => 'String',
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
    my %params = (
                  query_file => $self->query_file,
                  subject_files => \@subject_files,
                  psl_path => $self->psl_path,
                  blat_output_path => $self->blat_output_path
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
        $self->fatal_message('failed to run blat workflow.');
    }
    return 1;
}


1;

