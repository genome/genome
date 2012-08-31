package Genome::Model::Tools::PhredPhrap::PhdToFasta;

use strict;
use warnings;

use Genome;

require Cwd;
use Data::Dumper;

class Genome::Model::Tools::PhredPhrap::PhdToFasta {
    is => 'Command',
    has => [  
    phd_dir => {
        type => 'String',
        is_optional => 0,
        doc => 'Directory with PHDs',
    },
    phd_file => {
        type => 'String',
        is_optional => 0,
        doc => 'Files of PHDs to get FASTA and Quality',
    },
    fasta_file => {
        type => 'String',
        is_optional => 0,
        doc => 'FASTA output file.  Quality file will be named <FASTA>.qual',
    },
    _cwd => {
        type => 'String',
        is_optional => 1,
    },
    _error_file => {
        type => 'String',
        is_optional => 1,
    },
    ], 
};

sub help_brief {
    return 'Creates a FASTA and Qual file from PHDs';
}

sub help_detail {
    return 'Runs ewing_phd2fasta to create a FASTA and Qualty file from a list of PHDs in a file.  The FASTA file usually ends with a \'.fasta\' extension.  The Quality file will tahe the FASTA file name and add a \'.qual\' extension.';
}

sub qual_file { 
    my $self = shift;

    return sprintf('%s.qual', $self->fasta_file);
}

sub execute {
    my $self = shift;

    # Check phd dir
    $self->error_message(
        sprintf('Phd directory (%s) does not exist', $self->phd_dir)
    )
        and return unless -e $self->phd_dir;
    
    $self->error_message(
        sprintf('Phd directory (%s) is not a directory', $self->phd_dir)
    )
        and return unless -d $self->phd_dir;

    $self->error_message(
        sprintf('File of PHDs (%s) does not exist', $self->phd_file)
    )
        and return unless -e $self->phd_file;
    
    $self->error_message(
        sprintf('File of PHDs (%s) is empty', $self->phd_file)
    )
        and return unless -s $self->phd_file;

    $self->_error_file( $self->fasta_file . ".error" );
    unlink $self->_error_file if -e $self->_error_file;

    unlink $self->fasta_file if -e $self->fasta_file;
    unlink $self->qual_file if -e $self->qual_file;

    $self->_cwd( Cwd::getcwd() );
    chdir $self->phd_dir;

    my $command = sprintf(
        'ewing_phd2fasta -if %s -os %s -oq %s -of %s',
        $self->phd_file,
        $self->fasta_file,
        $self->qual_file,
        $self->_error_file,
    );

    my $rv = system $command;

    unlink $self->_error_file if -e $self->_error_file;
    
    chdir $self->_cwd;

    return ( $rv ) ? 0 : 1;
}

1;

=pod

=head1 Name

Genome::Model::Tools::PhredPhrap::PhdToFasta

=head1 Synopsis

Runs ewing_phd2fasta to create a FASTA and Qualty file from a list of PHDs in a file.  The FASTA file usually ends with a '.fasta' extension.  The Quality file will tahe the FASTA file name and add a '.qual' extension.

=head1 Usage

 use Genome;
 use Genome::Model::Tools::PhredPhrap::PhdToFasta

 my $phd2fnq =Genome::Model::Tools::PhredPhrap::PhdToFasta->create( # or execute
    phd_dir => 'PROJECT/phd_dir', # required
    phd_file => 'PROJECT/phds', # required 
    fasta_file => 'PROJECT/project.fasta, #required
 );

 $phd2fasta->execute;

=head1 Methods

=head2 execute

=over

=item I<Synopsis>   Runs the module

=item I<Arguments>  None

=item I<Returns>    True for success, False for failure

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
