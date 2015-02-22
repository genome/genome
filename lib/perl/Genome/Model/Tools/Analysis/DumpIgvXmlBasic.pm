package Genome::Model::Tools::Analysis::DumpIgvXmlBasic;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Cwd qw(abs_path);
use URI;
use XML::Writer;

class Genome::Model::Tools::Analysis::DumpIgvXmlBasic {
    is => 'Command::V2',
    has => [
        resource_files => {
            type => 'HASH',
            doc => 'A keyed list of files to display in the IGV session. Each key is a path to a file and the corresponding value is that label that should be displayed in the IGV session for that file.',
        },
        output_file => {
            type => 'String',
            is_optional => 0,
            doc => 'Output XML file',
        },
        reference_name => {
            type => 'String',
            is_optional => 0,
            doc => 'The name of the reference (in IGV) that the bams are aligned to.',
            example_values => ['b37', 'hg18', 'mm9']
        },
    ]
};

sub help_brief {
    "This helps create IGV session files for manual reviewers."
}

sub help_detail {
    <<'HELP';
Makes an IGV session for review. Supports as many tracks as you want. Panels and track layout are not explicitly specified.
HELP
}

sub execute {
    my $self=shift;

    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    my $xml_writer = new XML::Writer(
        OUTPUT => $output_fh,
        NEWLINES => 1,
        ENCODING => 'utf-8'
    );

    $xml_writer->xmlDecl('UTF-8', 'no');
    $xml_writer->startTag(
        'Session',
        genome => $self->reference_name,
        locus => '1:1-100',
        version => '3',
    );
    $xml_writer->startTag('Resources');
    while (my ($file, $label) = each %{$self->resource_files}) {
        $xml_writer->emptyTag(
            'Resource',
            path => resolve_file_path($file),
            relativePath => 'false',
            name => $label,
        );
    }
    $xml_writer->endTag('Resources');
    $xml_writer->endTag('Session');
    $xml_writer->end();
    $output_fh->close();

    return 1;
}

sub resolve_file_path {
    my $file = shift;
    my $uri = URI->new($file);
    my $path;
    if (defined($uri->scheme)) {
        if ($uri->scheme eq 'file') {
            $path = abs_path($uri->path);
        }
        else {
            $path = $uri->as_string;
        }
    }
    else {
        $path = abs_path($file);
    }
    return $path;
}

sub _is_hidden_in_docs {
    return 1;
}

1;
