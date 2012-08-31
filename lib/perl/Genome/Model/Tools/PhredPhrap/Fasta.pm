package Genome::Model::Tools::PhredPhrap::Fasta;

use strict;
use warnings;

use Genome;

require Cwd;
use Data::Dumper;
require File::Copy;

class Genome::Model::Tools::PhredPhrap::Fasta {
    is => 'Genome::Model::Tools::PhredPhrap',
    has => [
    fasta_file => {
        is => 'String', #file_r
        doc => "Fasta file.  If desired, have a quality file named '<FASTA>.qual'",
    },
    ],
};

sub help_brief {
    return 'Creates an assembly by running phrap on a FASTA (and Qual - <FASTA_FILE>.qual) file.';
}

sub help_detail {
    return help_brief();
}

sub memlog {
    my $self = shift;

    return sprintf('%s.memlog', $self->fasta_file);
}

sub out {
    my $self = shift;

    return sprintf('%s.phrap.out', $self->fasta_file);
}

sub execute {
    my $self = shift;

    for my $file (qw/ memlog out /) {
        unlink $self->$file if -e $self->$file;
    }

    my $cmd = sprintf('%s %s', $self->phrap_command_name, $self->fasta_file);
    my @properties = grep { $_->property_name ne 'version' } Genome::Model::Tools::PhredPhrap->get_class_object->direct_property_metas;
    for my $property ( @properties ) {
        my $property_name = $property->property_name;
        my $value = $self->$property_name;
        next unless defined $value;
        if ( $property->data_type eq 'Boolean' ) {
            next unless $value; # for 0
            $value = '';
        }
        $cmd .= sprintf(
            ' -%s %s', 
            $property_name, 
            $value,
        );
    }

    $cmd .= sprintf(' > %s 2> %s ', $self->out, $self->memlog);
    $self->status_message("Phrap: $cmd");
    my $phrap = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $phrap ) {
        $self->error_message("Phrap failed: $@");
        return;
    }

    return 1;
}

1;

=pod

=head1 Name

=head1 Synopsis

=head1 Methods

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This module is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
