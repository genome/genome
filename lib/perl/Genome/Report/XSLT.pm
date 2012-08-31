package Genome::Report::XSLT;
#:adukes check

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
use XML::LibXML;
use XML::LibXSLT;

class Genome::Report::XSLT {
};

sub transform_report {
    my ($class, %params) = @_;

    my $report = delete $params{report};
    unless ( $report ) {
        $class->error_message("Report is required to transform");
        return;
    }
               
    my $xslt_file = delete $params{xslt_file};
    my $validate = eval { Genome::Sys->validate_file_for_reading($xslt_file) };
    if (!$validate or $@) {
        $class->error_message("validate_file_for_reading on $xslt_file failed: $@");
        return;
    }

    my $xml = XML::LibXML->new();
    my $xml_string = $xml->parse_string( $report->xml_string );
    my $style_doc = eval{ $xml->parse_file( $xslt_file ) };
    if ( $@ ) {
        $class->error_message("Error parsing xslt file ($xslt_file):\n$@");
        return;
    }

    my $xslt = XML::LibXSLT->new();
    my $stylesheet = eval { $xslt->parse_stylesheet($style_doc) };
    if ( $@ ) {
        $class->error_message("Error parsing stylesheet for xslt file ($xslt_file):\n$@");
        return;
    }

    my $transformed_xml = $stylesheet->transform($xml_string);

    return {
        content => $stylesheet->output_string($transformed_xml),
        encoding => $stylesheet->output_encoding,
        media_type => $stylesheet->media_type,
        output_type => $class->_determine_output_type($stylesheet->media_type),
    };
}

sub _determine_output_type {
    my ($class, $media_type) = @_;

    if ( not defined $media_type ) {
        $class->error_message("No media type given to convert to output type");
        return '',
        #die instead??
    }
    if ( $media_type eq 'text/plain' ) {
        return 'txt';
    }
    elsif ( $media_type eq 'text/html' ) {
        return 'html';
    }
    elsif ( $media_type =~ /xml/ ) {
        return 'xml';
    }
    else {
        $class->error_message("Unknown media type to convert to output type: $media_type");
        return '';
        #die instead??
    }
}

1;

=pod

=head1 Name

Genome::Report::XSLT

=head1 Synopsis

=head1 Usage

 use Genome;
 
 # Get or generate a report...
 my $report = Genome::Report->create_report_from_directory(...);
 
 # Grab a xslt file
 my $xslt_file = ...;
 
 # Transform
 my $xslt = Genome::Report::XSLT->transform_report(
    report => $report, # required
    xslt_file => $xslt_file, #required
 );

 print "Content: ".$xslt->{content}."\n";

 ...
 
=head1 Public Methods

=head2 transform_report

 my $string = Genome::Report::XSLT->transform_report(report => $report, xslt_file => $xslt_file);

=over

=item I<Synopsis>   Takes a report and an xslt file (as a hash), and returns the transformed report as a string

=item I<Arguments>  report (Genome::Report), xslt_file (readable file)

=item I<Returns>    hashref with content, media_type, and encoding keys.

=back

=head1 Disclaimer

Copyright (C) 2009 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
