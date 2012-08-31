package Genome::Report::FromLegacy;
#:adukes check

use strict;
use warnings;

use Genome;

require File::Basename;
use Storable;

class Genome::Report::FromLegacy {
    is => 'Genome::Report::Generator',
    has => [
    name => {
        is => 'Text',
        doc => 'Name to give the report.  Will usually have a default/calculated value',
    },
    description => {
        is => 'Text',
        doc => 'Description to give the report.  Will usually have a default/calculated value',
    },
    ],
};

sub create {
    my ($class, %params) = @_;

    my $properties_file = delete $params{properties_file};
    unless ( Genome::Sys->validate_file_for_reading($properties_file) ) {
        $class->error_message("Properties file not given or does not exist");
        return;
    }

    my $properties = Storable::retrieve($properties_file);
    unless ( $properties ) {
        $class->error_message(
            sprintf(
                "Can't get properties from legacy file (%s)\n%s",
                $properties_file,
                $!,
            )
        );
        return;
    }
    $properties->{generation_params} ||= {};

    my $dir = File::Basename::dirname($properties_file);
    $dir =~ s#/$##;
    $params{name} = Genome::Report->directory_to_name($dir);
    $params{description} = delete $properties->{description};
    
    my $self = $class->SUPER::create(%params)
        or return;

    $self->{_properties} = $properties;
    
    return $self;
}
    
sub _add_to_report_xml { 
    return 1; 
}

sub _properties {
    return $_[0]->{_properties};
}

sub _get_date {
    return $_[0]->_properties->{date} || $_[0]->SUPER::_get_date;
}

sub _get_generator {
    return $_[0]->_properties->{generator};
}

sub _get_params_for_generation {
    return %{$_[0]->_properties->{generation_params}};
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

