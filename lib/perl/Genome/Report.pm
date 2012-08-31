package Genome::Report;
#:adukes per eddie's suggestion, this could use a method to modify and resave report data

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
require File::Basename;
use Storable;
require XML::LibXML;

class Genome::Report {
    is => 'UR::Object',
    has => [
    xml => { is => 'XML::LibXML::Document', doc => 'XML document object', },
    parent_directory => {
        is => 'Text',
        is_optional => 1,
        doc => 'Parent directory where report lives.',
    },
    name => {
        calculate => q| return $self->_get_report_attribute('name'); |,
        doc => 'Report name',
    },
    description => {
        calculate => q| return $self->_get_report_attribute('description'); |,
        doc => 'Report description',
    },
    date => {
        calculate => q| return $self->_get_report_attribute('date'); |,
        doc => 'Date the report was generated',
    },
    generator => {
        calculate => q| return $self->_get_report_attribute('generator'); |,
        doc => 'The class of the generator that created this report',
    },
    ],
};

#sub xml { my ($self, $xml) = @_; $self->{_xml} = $xml if $xml; return $self->{_xml}; }
#< XML Attributes >#
sub xml_string {
    return $_[0]->xml->toString(1);
}

sub xml_root_element {
    my $self = shift;
    my $root = $self->xml->documentElement; 
    confess "Can't get root element from XML Document\n" unless $root;
    return $root;
}

sub _get_report_attribute { 
    my ($self, $attr) = @_;
    my ($attribute_element) = $self->xml->findnodes('report/report-meta/'.$attr);
    # confessing here cuz this is a private method and should always return something
    confess "Can't get attribute element for $attr from XML Document\n" unless $attribute_element;
    return $attribute_element->to_literal;
}

#< Datasets >#
sub get_dataset_nodes { 
    return $_[0]->xml->findnodes('report/datasets/*');
}

sub get_dataset_nodes_for_name { 
    my ($self, $name) = @_;

    confess "No name given to get dataset node by name." unless $name;

    return $_[0]->xml->findnodes("report/datasets[1]/$name");
}

sub get_dataset {
    my ($self, $name) = @_;

    my ($node) = $self->get_dataset_nodes_for_name($name)
        or return;

    return Genome::Report::Dataset->create_from_xml_element($node);
}

sub get_datasets {
    my ($self, $name) = @_;

    return map { 
        Genome::Report::Dataset->create_from_xml_element($_) 
    } $self->get_dataset_nodes;
}

sub get_dataset_names { 
    my $self = shift;

    my @datasets = $self->get_dataset_nodes
        or return;
    
    my @names;
    for my $dataset ( @datasets ) {
        push @names, $dataset->nodeName;
    }
    
    return @names;
}

sub get_dataset_by_name_as_xml {
    my ($self, $name) = @_;

    my ($dataset) = $self->get_dataset_nodes_for_name($name)
        or return;

    return $dataset->toString;
}

sub get_datasets_by_name_as_separated_value_string {
    my ($self, $name, $separator) = @_;

    my $dataset = $self->get_dataset($name)
        or return;

    return $dataset->to_separated_value_string(separator => $separator);
}

#<>#
sub generator_params {
    my $self = shift;

    my %params;
    for my $attr_element ( $self->xml->findnodes('report/report-meta/generator-params/*') ) {
        my $attr = $attr_element->nodeName;
        $attr =~ s#\-#\_#g;
        push @{$params{$attr}}, $attr_element->to_literal;
    }
    #print Dumper(\%params);

    return \%params;
}

#< DATA - BACKWARD COMPATIBILITY - THIS WILL BE REMOVED! >#
sub data {
    my ($self, $data) = @_;

    $self->{_data} = $data if defined $data;

    return $self->{_data};
}

#< Get/Create >#
sub get { # no real 'get'...since storage is not synced
    confess "Cannot conventionally 'get' a report.  To get a stored report, use 'create_report_from_directory' or 'get_reports_in_parent_directory'\n";
}

sub create {
    my ($class, %params) = @_;

    my $data = delete $params{data};

    my $xml = delete $params{xml};
    my $self = $class->SUPER::create(%params)
        or return;
    $self->xml($xml);

    unless ( $self->xml ) {
        $self->error_message("Need XML document object to create report");
        $self->delete;
        return;
    }

    return $self;
}

sub create_report_from_directory {
    my ($class, $directory) = @_;

    Genome::Sys->validate_directory_for_read_access($directory)
        or return;

    $class->_convert_properties_file_to_xml($directory);
    
    my %params;
    $params{xml} = $class->_retrieve_xml($directory)
        or return;

    $params{parent_directory} = File::Basename::dirname($directory);
    
    return $class->create(%params);
}

sub create_reports_from_parent_directory {
    my ($class, $parent_directory) = @_;

    Genome::Sys->validate_directory_for_read_access($parent_directory)
        or return;

    my @reports;
    for my $directory ( glob($parent_directory.'/*') ) {
        next unless -d $directory;
        my $report;
        eval {
            $report = $class->create_report_from_directory($directory)
        };
        next unless $report;
        push @reports, $report;
    }

    return @reports;
}

#< File and Directory Naming >#
sub directory {
    my $self = shift;

    # Validate parent dir
    my $parent_directory = $self->parent_directory;
    unless ( $parent_directory ) {
        $self->error_message("No parent directory set for report: ".$self->name);
        return;
    }
    my $validate = eval {
        Genome::Sys->validate_directory_for_write_access( 
            $parent_directory
        ); 
    };
    if (!$validate or $@) {
        $self->error_message("validate_directory_for_write_access on $parent_directory failed: $@");
        return;
    }

    # Get dir, create
    my $directory = $parent_directory.'/'.join('_', split(' ', $self->name));
    Genome::Sys->create_directory($directory)
        or return;

    return $directory;
}

sub _xml_file {
    my ($class, $directory) = @_;

    confess "Need directory to get xml file name\n" unless $directory;
    
    return $directory.'/report.xml';
}

sub name_to_subdirectory {
    return join('_', split(' ', $_[1]));
}

sub directory_to_name {
    my ($class, $directory) = @_;

    $directory =~ s#/$##;
    my $basename = File::Basename::basename($directory);
    return join(' ', split('_', $basename));
}

#< Save/Retrieve >#
sub save {
    my ($self, $parent_directory, $overwrite) = @_;

    unless ( $self->xml ) {
        $self->error_message('No XML was found to save');
        return;
    }

    $self->parent_directory($parent_directory);
    my $directory = $self->directory
        or return;

    # File
    my $file = $self->_xml_file($directory)
        or return;
    if ( $overwrite ) {
        unlink $file if -e $file;
    }
    my $validate = eval { Genome::Sys->validate_file_for_writing($file)};
    if (!$validate or $@) {
        $self->error_message("validate_file_for_writing on $file failed: $@");
        return;
    }

    # Save it
    unless ( $self->xml->toFile($file, 1) ) {
        $self->error_message("Could not save report xml to file ($file): $!");
        return;
    }

    # DATA - BACKWARD COMPATIBILITY - THIS WILL BE REMOVED!
    my $data = $self->data;
    if ( $data ) {
        for my $type (qw/ csv html txt /) {
            next unless exists $data->{$type};
            my $file = $directory.'/report.'.$type;
            unlink $file if -e $file;
            my $fh = Genome::Sys->open_file_for_writing($file)
                or confess;
            $fh->print( $data->{$type} );
            $fh->close;
        }
    }

    return 1;
}

sub _retrieve_xml {
    my ($class, $directory) = @_;

    my $file = $class->_xml_file($directory);
    Genome::Sys->validate_file_for_reading($file)
        or return;

    my $libxml = XML::LibXML->new();
    my $xml;
    eval {
        $xml = $libxml->parse_file($file);
    };
    unless ( $xml ) {
        $class->error_message("Could not parse report XML from file ($file): $@");
        return;
    }

    return $xml;
}

sub _convert_properties_file_to_xml { # for legacy reports
    my ($class, $directory) = @_;
    
    my $properties_file = $directory.'/properties.stor';
    return 1 unless -s $properties_file; # nothing to do
    
    my $xml_file = $class->_xml_file($directory); 
    if ( -e $xml_file ) { # has both, just remove props file
        unlink $properties_file;
        return 1;
    }
    
    # convert
    my $generator = Genome::Report::FromLegacy->create(properties_file => $properties_file);
    unless ( $generator ) {
        $class->error_message("Can't create report from legacy generator");
        return;
    }
    my $report = $generator->generate_report;
    unless ( $report ) {
        $class->error_message("Can't generate report to convert legacy report");
        return;
    }

    # save
    $directory =~ s#/$##; # just in case, file::basename don't like trailers
    my $subdir = File::Basename::dirname($directory);
    $report->save($subdir, 1);

    # remove props file
    unlink $properties_file;

    return 1;
}

# Old data methods
sub get_brief_output {
    return $_[0]->description;
    # FYI this was generating if file didn't exist
}

sub get_detail_output {
    die;
    return $_[0]->get_as_html;
    # FYI this was generating if file didn't exist
}

1;

=pod

=head1 Name

Genome::Report

=head1 Synopsis

A generic report object that stores on the filesystem.

=head1 Usage

 my $report = Genome::Report->create(
    name => 'Happy', # required
    data => { # required
        generator => 'Genome::Report::Happy', # required
        generator_params => { # required
            happy_level => 9,
        },
        date => 'today',
        description => 'This is a happy report', # required
        html => '<html></html>', # html data
        csv => ..., # csv data
        xml => ..., # xml data
    },
 );

 $report->save('some_directory')
    or die;

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
