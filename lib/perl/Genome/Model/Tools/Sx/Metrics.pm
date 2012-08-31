package Genome::Model::Tools::Sx::Metrics;

use strict;
use warnings;

use Genome;

use XML::LibXML;
use XML::LibXSLT;

class Genome::Model::Tools::Sx::Metrics {
    is_abstract => 1,
};

sub metric_names { Carp::confess('Define metric names in subclass!'); }
sub calculate { return 1; }

sub get_metric {
    my ($self, $name) = @_;

    return $self->$name;
}

sub set_metric {
    my ($self, $name, $value) = @_;

    return $self->$name($value);
}

sub to_xml {
    my $self = shift;
    $self->calculate;
    my $view = UR::Object::View::Default::Xml->create(subject => $self);
    return $view->content;
}

sub transform_xml {
    my ($self, $xslt_file) = @_;

    Carp::confess('No xslt file given to transform xml!') if not $xslt_file;

    my $xml = XML::LibXML->new();
    my $xml_string = $xml->parse_string( $self->to_xml );
    my $style_doc = eval{ $xml->parse_file( $xslt_file ) };
    if ( $@ ) {
        Carp::confess("Error parsing xslt file ($xslt_file):\n$@");
    }

    my $xslt = XML::LibXSLT->new();
    my $stylesheet = eval { $xslt->parse_stylesheet($style_doc) };
    if ( $@ ) {
        Carp::confess("Error parsing stylesheet for xslt file ($xslt_file):\n$@");
    }

    my $transformed_xml = $stylesheet->transform($xml_string);
    return $stylesheet->output_string($transformed_xml);
}

sub transform_xml_to {
    my ($self, $type) = @_;
    my $xslt_file = $self->xslt_file_for($type);
    return $self->transform_xml($xslt_file);
}

sub xslt_file_for {
    my ($self, $type) = @_;

    Carp::confess('No type to get xslt file!') if not $type;
    
    my $genome_dir = Genome->get_base_directory_name;
    #my $inc_dir = substr($genome_dir, 0, -7); # rm Genome
    my $module = $self->class;
    $module =~ s/Genome:://;
    $module =~ s#::#/#g;
    my $xslt_file = sprintf(
        '%s/%s.%s.xsl',
        $genome_dir,
        $module,
        'txt',#$type
    );
    if ( not -s $xslt_file ) {
        Carp::confess("No xslt file ($xslt_file) for type ($type)");
    }

    return $xslt_file;
}

sub to_text {
    my $self = shift;
    $self->calculate;
    my $view = UR::Object::View::Default::Text->create(subject => $self);
    return $view->content;
}

sub from_file {
    my ($class, $file) = @_;

    if ( not $file ) {
        $class->error_message('No file given to create metrics from file');
        return;
    }

    if ( not -s $file ) {
        $class->error_message('Failed to read metrics from file. File ('.$file.') does not exist.');
        return;
    }

    my $fh = eval{ Genome::Sys->open_file_for_reading($file); };
    if ( not $fh ) {
        $class->error_message("Failed to open file ($file)");
        return;
    }

    my $line = $fh->getline;
    my $metrics_class = 'Genome::Model::Tools::Sx::Metrics::Basic';
    my $self;
    if ( $line =~ /Sx Metrics/ ) { # class name
        chomp $line;
        my ($id) = $line =~ s/ '(.+)'$//; # id w/ spaces is in quotes, grab it
        my @tokens = split(/\s+/, $line);
        $id = pop @tokens if not $id; # if id is not in quotes, pop it off the end
        $metrics_class = join('::', @tokens);
        $metrics_class = 'Genome::Model::Tools::'.$metrics_class if not $metrics_class =~ /^Genome::Model::Tools::/;
        $self = $metrics_class->get($id); # it may be in the system already
        $self = $metrics_class->create(id => $id) if not $self; # create it if get does not work
    }
    else {
        $fh->seek(0, 0); # legacy metrics file does not have a class; reset fh to begining
        $self = $metrics_class->create() if not $self; # new...create!
    }

    if ( not $self ) {
        $class->error_message("Failed to create metrics ($metrics_class)!");
        return;
    }

    while ( $line = $fh->getline ) {
        chomp $line;
        my ($key, $val);
        if ( $line =~ /=/ ) {
            ($key, $val) = split('=', $line, 2) if not $key; # legacy
        }
        else {
            ($key, $val) = split(': ', $line, 2);
        }
        $key =~ s/^\s+//;
        $self->$key($val);
    }
    $fh->close;

    return $self;
}

sub to_file {
    my ($self, $file) = @_;

    if ( not $file ) {
        $self->error_message('No file given to create metrics from file');
        return;
    }

    unlink $file;
    my $fh = eval{ Genome::Sys->open_file_for_writing($file); };
    if ( not $fh ) {
        $self->error_message("Failed to open file ($file)");
        return;
    }
    $fh->print($self->to_text);
    $fh->close;

    return 1;
}

1;

