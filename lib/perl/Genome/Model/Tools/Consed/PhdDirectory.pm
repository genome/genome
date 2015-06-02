package Genome::Model::Tools::Consed::PhdDirectory;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Consed::PhdDirectory {
    has => [
        directory => {
            is => 'Text',
        },
    ],
};

#- GENERAL -#
sub base_names {
    my $self = shift;

    my %names;
    for my $phd ( $self->phd_names ) {
        $phd =~ s/\.phd\.\d+$//;
        $names{$phd} = 1;
    }

    return keys %names;
}

sub phd_names {
    my $self = shift;

    my $directory = $self->directory;
    my @phd_files = glob("$directory/*");

    return grep { s#$directory/##; m#\.phd\.# } @phd_files;
}

sub phd_file {
    my ($self, $name) = @_;

    return sprintf('%s/%s', $self->directory, $name);
}

sub phd_file_name {
    return phd_file(@_);
}

sub phd_exists
{
	my ($self, $name) = @_;

    my $file = $self->phd_file_name($name);

	return unless -s $file;

    return $file;
}

sub get_phd
{
	my ($self, $phd_name) = @_;
	
    my $file = $self->phd_exists($phd_name)
        or return;

    my $phd = Genome::Model::Tools::Consed::PhdReader->read($file);
    if ( not $phd ) {
        $self->error_message('No phd frpm reader for file: '.$file);
        return;
    }

    return $phd;
}

#- LATEST PHD -#
sub latest_phd_iteration
{
    my ($self, $name) = @_;

    my $latest_iteration = 0;
    foreach my $phd_name ( glob( sprintf('%s/%s.phd.*', $self->directory, $name) ) )
    {
        my ($iteration) = $phd_name =~ /\.phd\.(\d+)$/;
        $latest_iteration = $iteration if $iteration > $latest_iteration;
    }

    return $latest_iteration;
}

sub latest_phd_name
{
    my ($self, $name) = @_;

    my $latest_iteration = $self->latest_phd_iteration($name);
    
    return if $latest_iteration == 0;
    
	return sprintf('%s.phd.%s', $name, $latest_iteration);
}

sub latest_phd_file
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;

    return $self->phd_file_name($phd_name);
}

sub latest_phd
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;
    
    return $self->get_phd($phd_name);
}

sub get_latest_phd
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;
    
    return $self->get_phd($phd_name);
}

#- NEXT PHD -#
sub next_phd_name
{
    my ($self, $name) = @_;

    return sprintf('%s.phd.%d', $name, $self->latest_phd_iteration($name) + 1); 
}

sub next_phd_file
{
    my ($self, $name) = @_;

    return $self->phd_file_name( $self->next_phd_name($name) );
}

1;		

