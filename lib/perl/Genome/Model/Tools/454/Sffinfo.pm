package Genome::Model::Tools::454::Sffinfo;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::454::Sffinfo {
    is => ['Genome::Model::Tools::454'],
    has => [
            sff_file => {
                         doc => 'The sff file to operate',
                         is => 'string',
                     },
            output_file => {
                            is => 'string',
                            doc => 'The output file path',
                        },
            params => {
                       is => 'string',
                       doc => 'The params to pass to sffinfo',
                   },
        ],
};

sub help_brief {
    "convert sff file to fasta file"
}

sub help_detail {
    return <<EOS
convert sff file to fasta file
see sffinfo usage for valid params
EOS
}
sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('This genome-model tool '. $self->command_name .' will only run on 64-bit');
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $cmd = $self->bin_path .'/sffinfo '. $self->params .' '. $self->sff_file .' > '. $self->output_file;
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero exit code '$rv' returned by sffinfo");
        return;
    }
    return 1;
}

1;


