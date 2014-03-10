package Genome::Model::Tools::PhredPhrap::AbiToScf;

use strict;
use warnings;
  
use Genome;

use Data::Dumper;
require File::Copy;
require IO::Dir;

class Genome::Model::Tools::PhredPhrap::AbiToScf {
    is => 'Command',
    has => [
    abi_dir => {
        is => 'String', #dir_r,
        is_optional => 0,
        doc => 'Directory of ABIs.',
    },
    abi_file => {
        is => 'String', #dir_r,
        is_optional => 1,
        doc => 'File of ABIs to convert to SCFs.',
    },
    chromat_dir => { 
        is => 'String', #dir_rw
        is_optional => 0,
        doc => 'Directory to put SCFs',
    },
    scf_ext => {
        is => 'String',
        is_optional => 1,
        default => 'scf',
        doc => 'Extention to replace ab[i1] with in new SCF file. Default is \'scf\'.  Use NONE for no extention',
    },
    ],
};

sub help_brief {
    return 'Runs phred-gasp to convert ABIs to SCFs.';
}

sub help_detail {
    return help_brief();
}

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_)
        or return;

    if ($self->scf_ext eq 'NONE') {
      $self->scf_ext('');
    }

    for my $dir (qw/ chromat_dir abi_dir /) {
        unless ( -e $self->$dir ) { 
            $self->error_message( sprintf('%s (%s) does not exist', $dir, $self->$dir) );
            return;
        }

        unless ( -d $self->$dir ) { 
            $self->error_message( sprintf('%s (%s) is not a directory', $dir, $self->$dir) );
            return;
        }
    }

    if ( $self->chromat_dir eq $self->abi_dir ) { 
        $self->error_message("Cannot have same directory for abis and SCFs");
        return;
    }

    return $self;
}

sub execute {
    my $self = shift;
    
    my %phred_gasp_params = (
        cv => 3,
        cp => 2,
        cd => $self->chromat_dir,
    );

    if ( $self->abi_file ) {
        $phred_gasp_params{'if'} = $self->abi_file;
    }
    else {
        $phred_gasp_params{id} = $self->abi_dir;
    }
    
    my $cmd = sprintf('phred-gasp %s -c', join(' ', map { sprintf('-%s %s', $_, $phred_gasp_params{$_}) } keys %phred_gasp_params));
    $self->debug_message("RUNNING:\n $cmd");
    if ( system($cmd) ) {
        $self->error_message("Failed to run phred-gasp: $@");
        return;
    }
    
    my $io_dir = IO::Dir->new( $self->chromat_dir );
    unless ( $io_dir ) {
        $self->error_message( sprintf('Can\'t open chromat_dir (%s): %s', $self->chromat_dir) );
        return;
    }

    my $new_ext = ( defined $self->scf_ext and $self->scf_ext ne '' )
    ? sprintf('.%s', $self->scf_ext)
    : '';

    while ( my $file = $io_dir->read ) {
        chomp $file;
        my $abi_file = sprintf('%s/%s', $self->chromat_dir, $file);
        next unless $file =~ s/\.ab[i1]$/$new_ext/;
        my $scf_file = sprintf('%s/%s', $self->chromat_dir, $file);
        $self->error_message("Can't copy $abi_file to $scf_file\: $!") unless File::Copy::move($abi_file, $scf_file);
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
