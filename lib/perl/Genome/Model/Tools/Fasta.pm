package Genome::Model::Tools::Fasta;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
require File::Copy;

class Genome::Model::Tools::Fasta {
    is => 'Command',
    has_input => [
        fasta_file => {
            type => 'String',
            is_optional => 0,
            doc => 'FASTA file. Quality file (if appropriate) will be named <fasta_file>\'.qual\'',
        },
    ],
};

sub help_brief {
    return "Tools for working with FASTA and Qual files"
}

sub help_detail {
    return help_brief();
}

sub create { 
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;
    
    unless ( defined $self->fasta_file ) {
        $self->error_message("A fasta file is required");
        return;
    }

    my $validate = eval { Genome::Sys->validate_file_for_reading( $self->fasta_file ); };
    if (!$validate or $@) {
        $self->error_message("validate_file_for_reading failed: $@");
        return;
    }
    
    my ($basename, $directory, $suffix) = File::Basename::fileparse($self->fasta_file, '.fasta', '.fas', '.fa', '.fna');
    unless ( $suffix ) {
        $self->error_message( sprintf('FASTA file (%s) needs to have a ".fasta", ".fas" or ".fa" suffix.', $self->fasta_file) );
        return;
    }

    $self->{_cwd} = Cwd::getcwd();
    $self->{_fasta_directory} = $directory;
    $self->{_fasta_basename} = $basename;
    $self->{_fasta_suffix} = $suffix; # Remember it has the '.' in it!

    return $self;
}

sub DESTROY {
    my $self = shift;

    $self->chdir_cwd if defined $self->{_cwd};
    
    return 1;
}

#< Dirs >#
sub cwd {
    return $_[0]->{_cwd};
}

sub chdir_cwd {
    my $self = shift;

    unless (  chdir $self->cwd ) {
        $self->error_message( sprintf('Can\'t access cwd (%s)', $self->cwd) );
        return;
    }

    return 1;
}

sub chdir_fasta_directory {
    my $self = shift;

    unless (  chdir $self->fasta_directory ) {
        $self->error_message( sprintf('Can\'t access fasta directory (%s)', $self->fasta_directroy) );
        return;
    }

    return 1;
}

#< FASTA Fileparse #>
sub fasta_basename {
    return $_[0]->{_fasta_basename};
}

sub fasta_directory {
    return $_[0]->{_fasta_directory};
}

sub fasta_suffix {
    return $_[0]->{_fasta_suffix};
}

sub fasta_base {
    return sprintf('%s%s', $_[0]->fasta_basename, $_[0]->fasta_suffix);
}

#< Qual file >#
sub qual_base {
    return sprintf('%s.qual', $_[0]->fasta_base);
}

sub qual_file {
    return sprintf('%s%s', $_[0]->fasta_directory, $_[0]->qual_base);
}

sub have_qual_file {
    return -s $_[0]->qual_file;
}

#< New file names >#
sub fasta_base_with_new_suffix { 
    my ($self, $suffix) = @_;

    return sprintf('%s.%s%s', $self->fasta_basename, $suffix, $self->fasta_suffix);
}

sub fasta_file_with_new_suffix { 
    my ($self, $suffix) = @_;

    return sprintf('%s%s', $self->fasta_directory, $self->fasta_base_with_new_suffix($suffix));
}

sub qual_base_with_new_suffix {
    my ($self, $suffix) = @_;

    return sprintf('%s.qual', $self->fasta_base_with_new_suffix($suffix));
}

sub qual_file_with_new_suffix {
    my ($self, $suffix) = @_;

    return sprintf('%s%s', $self->fasta_directory, $self->qual_base_with_new_suffix($suffix));
}

#< Back Up >#
sub default_back_up_suffix {
    return 'bak';
}

sub fasta_back_up_file {
    my ($self, $ext) = @_;

    return $self->fasta_file_with_new_suffix( defined $ext ? $ext : $self->default_back_up_suffix );
}

sub qual_back_up_file {
    my ($self, $ext) = @_;

    return $self->qual_file_with_new_suffix( defined $ext ? $ext : $self->default_back_up_suffix );
}

sub back_up_fasta_and_qual_files {
    my ($self, $ext) = @_;

    $ext = $self->default_back_up_suffix unless defined $ext;

    my $fasta_bak = $self->back_up_fasta_file($ext)
        or return;

    my $qual_bak = $self->back_up_qual_file($ext)
        or return;

    return ( $fasta_bak, $qual_bak );
}

sub back_up_fasta_file {
    my ($self, $ext) = @_;

    my $fasta_bak = $self->fasta_back_up_file($ext);
    unlink $fasta_bak if -e $fasta_bak;

    unless ( File::Copy::copy($self->fasta_file, $fasta_bak) ) {
        $self->error_message( sprintf('Can\'t copy %s to %s', $self->fasta_file, $fasta_bak) );
        return;
    }

    return $fasta_bak;
}

sub back_up_qual_file {
    my ($self, $ext) = @_;

    my $qual_bak = $self->qual_back_up_file($ext);
    unlink $qual_bak if -e $qual_bak;

    unless ( File::Copy::copy($self->qual_file, $qual_bak) ) {
        $self->error_message( sprintf('Can\'t copy %s to %s', $self->qual_file, $qual_bak) );
        return;
    }

    return $qual_bak;
}

#< Bio::SeqIO stuff >#
sub get_fasta_reader {
    return _get_bioseq_reader(@_, 'fasta');
}

sub get_qual_reader {
    return _get_bioseq_reader(@_, 'qual');
}

sub _get_bioseq_reader {
    Genome::Sys->validate_file_for_reading($_[1]) 
        or return;
    return _get_bioseq(@_, '<');
}

sub get_fasta_writer {
    Genome::Sys->validate_file_for_writing($_[1]) 
        or return;
    return _get_bioseq_writer(@_, 'fasta');
}

sub get_qual_writer {
    Genome::Sys->validate_file_for_writing($_[1]) 
        or return;
    return _get_bioseq_writer(@_, 'qual');
}

sub _get_bioseq_writer {
    return _get_bioseq(@_, '>');
}

sub _get_bioseq {
    my ($self, $file, $format, $rw) = @_;

    return Bio::SeqIO->new(
        '-file' => $rw.' '.$file,
        '-format' => $format,
    );
}


1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

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

