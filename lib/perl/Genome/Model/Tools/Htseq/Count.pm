package Genome::Model::Tools::Htseq::Count;
use strict;
use warnings;
use Genome;

use File::Spec;

class Genome::Model::Tools::Htseq::Count {
    is => 'Genome::Command::DelegatesToResult',
    has_param => [
        app_version => {
            is => 'SoftwareVersion',
            default_value => '0.5.4p1',
            valid_values => [ Genome::Sys->sw_versions('htseq','htseq-count') ],
            is_param => 1,
            doc => 'the version of htseq-count to use',
        },
        result_version => {
            is => 'Integer',
            default_value => 2,
            valid_values => [ 1, 2 ],
            doc => 'the version of results, which may iterate as this logic iterates',
        },
        mode => {
            is => 'Text',
            default_value => 'intersection-strict',
            valid_values => [ "intersection-strict" ],
            doc => 'mode',
        },
        minaqual => {
            is => 'Integer',
            default_value => 1,
            doc => 'minaqual',
        },
        whitelist_alignments_flags => {
            is => 'Text',
            is_optional => 1,
            is_param => 1,
            doc => 'require alignments to match the specified flags (-f): 0x0002 limits to only properly-paired alignments',
        },
        blacklist_alignments_flags => {
            is => 'Text',
            default_value => '0x0104',
            is_optional => 1,
            is_param => 1,
            doc => 'exclude alignments which match the specified flags (-F): 0x0104 excludes non-primary alignments and unaligned reads',
        },
        limit => {
            is => 'Number',
            is_optional => 1,
            is_param => 1,
            doc => 'limit the number of alignments to the first N (for testing)',
        },
        #sort_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION [-p1 -p2]'
        #},
        #bam_view_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION [-p1 -p2]'
        #},
    ],
    has_input => [
        alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            where => [ 'instrument_data.sample.extraction_type in' => [ "rna", "cdna", "total rna" ] ],
            example_values => [ { "instrument_data.sample.individual.common_name like" => "HCC%" } ],
            is_many => 1,
            doc => 'alignment results, typically from an RNA aligner',
        },
    ],
    has_optional => [
        output_dir => {
            is => 'FilesystemPath',
            doc => 'Location to symlink the result output_dir',
        },
    ],
    doc => 'generate htseq results for an (annotation-based) alignment result',
};

sub result_class {
    return __PACKAGE__ . '::Result';
}

sub help_synopsis {
    return <<EOS

gmt htseq count --alignment-results "instrument_data.id=2890686892" --app-version 0.5.4p1

gmt htseq count --alignment-results "instrument_data.sample.name='H_MU-752713-1209062'" --app-version 0.5.4p1

gmt htseq count --alignment-results "instrument_data.sample.individual.common_name like 'HCC%'" --app-version 0.5.4p1

# skip any data sets flagged as test data
gmt htseq count --alignment-results "instrument_data.sample.individual.common_name like 'HCC%' and test_name is null"
EOS
}

sub help_detail {
    return <<EOS
This tool runs "htseq-count" from the HTSeq package, developed by Simon Anders at EMBL Heidelberg (Genome Biology Unit).
http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
EOS
}

sub _additional_help_sections {
    return (
       "INPUTS",
       <<EOS,
It operates on "alignment results" from RNA/cDNA instrument data, such as those produced by tophat or rna-star.  These
results have associated annotation used during alignment, and that annotation is fed into htseq.

When run on more than one alignment result, each is processed individually, and the results are also aggregated into an additinal data product which is returned.
EOS
       "OUTPUTS",
       <<EOS,
The output is a list if "hit counts" per gene, and a second list of hit counts per transcript.
EOS
        "NOTE",
        <<EOS,
This tool saves software results with each execution, and shortcuts on subsequent runs to avoid duplicating effort.
EOS
  );
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
    # expect to return POD
}

sub _doc_authors {
    return <<EOS
 Scott Smith
 Malachi Griffith, Ph.D.
 Obi Griffith, Ph.D.
EOS
}

sub _doc_copyright_years {
    (2013);
}

sub _doc_license {
    my $self = shift;
    my (@y) = $self->_doc_copyright_years;  
    my $range;
    if (@y == 1) {
        $range = "$y[0]";
    }
    elsif (@y > 1) {
        $range = "$y[0]-$y[-1]";
    }
    return <<EOS
Copyright (C) $range Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the 
associated LICENSE file in this distribution.
EOS
}

sub _doc_credits {
    return ('','Simon Anders at EMBL Heidelberg (Genome Biology Unit) is the author of HTSeq.');
}

sub _doc_see_also {
    return <<EOS
B<Genome::Model::RnaSeq>(3), B<Genome::InstrumentData::AlignmentResult::PerLaneTophat>(3)
B<genome-model-rnaseq>(1), B<genome-instrument-data-align-per-lane-tophat>(1)
EOS
}

sub post_get_or_create {
    my $self = shift;
    my $result = $self->output_result;

    if($self->output_dir) {
        my ($vol, $dir, $file) = File::Spec->splitpath($self->output_dir);
        Genome::Sys->create_directory($dir);
        Genome::Sys->create_symlink($result->output_dir, $self->output_dir);
    }

    return 1;
}

1;

