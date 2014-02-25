package Genome::Model::Tools::Park::Base;

use strict;
use warnings;

use Genome;
use File::Basename qw(dirname);
use Carp qw(confess);

class Genome::Model::Tools::Park::Base {
    is => 'UR::Object',
    is_abstract => 1,
};

sub _template_path {
    die 'abstract',
}

sub _generate_inputs_file {
    my $self = shift;

    my $template_fh = Genome::Sys->open_file_for_reading($self->_template_path);
    my ($inputs_fh, $inputs_filename) = Genome::Sys->create_temp_file();

    my $line_number = 0;
    while (my $line = $template_fh->getline) {
        $line_number++;
        my ($workflow_input_name, $accessor) = _parse_line($line);

        unless ($self->can($accessor)) {
            confess(sprintf("Class (%s) has no accessor (%s) that the template (%s) expected on line (%d)",
                    ref $self, $accessor, $self->_template_path, $line_number));
        }
        if (!defined($self->$accessor)) {
            confess(sprintf("Attribute (%s) is undefined, which is disallowed in an inputs file.", $accessor));
        }
        my $value = [$self->$accessor] if ref($self->$accessor) ne 'ARRAY';

        my $string = join("\t", @$value);
        printf $inputs_fh "%s\t%s\n", $workflow_input_name, $string;
    }
    return $inputs_filename;
}

sub _parse_line {
    my $line = shift;

    chomp($line);
    my @parts = split(/\t/, $line);
    my $workflow_input_name = shift(@parts);
    my $accessor = shift(@parts);

    unless (defined($accessor)) {
        $accessor = (split(/\./, $workflow_input_name))[-1];
    }
    return ($workflow_input_name, $accessor);
}

1;
