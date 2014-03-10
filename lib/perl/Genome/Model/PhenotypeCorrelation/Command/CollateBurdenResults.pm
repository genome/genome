package Genome::Model::PhenotypeCorrelation::Command::CollateBurdenResults;

use Genome;
use Carp qw(confess);
use JSON;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::CollateBurdenResults {
    is => "Command::V2",
    has_input => [
        results_directory => {
            is => "FilePath",
            doc => "The directory containing the burden test results",
        },
        output_file => {
            is => "FilePath",
            doc => "Where to write the collated results",
        },
        job_params => {
            is => "Text",
            is_many => 1,
            doc => "list of JSON encoded hashrefs containing gene, phenotype, analysis_data_type, covariates",
        },
        dummy => {
            is_many => 1,
            is => "Text",
            doc => "Used to force workflow to block until burden test is done",
            is_optional => 1,
        }
    ],
    has_transient_optional => [
        _header_source_file => {
            is => "Text",
        },
        _header_fields => {
            is => "ARRAY",
        },
    ],
};

sub _parse_result_file {
    my ($self, $result_file) = @_;

    my $fh = Genome::Sys->open_file_for_reading($result_file);
    my $header_line = $fh->getline;
    chomp $header_line;
    my @header_fields = split(",", $header_line);
    my $stored_header_fields = $self->_header_fields;
    if ($stored_header_fields) {
        if (scalar @header_fields != scalar @$stored_header_fields) {
            my $header_src = $self->_header_source_file;
            confess $self->error_message(
                "Header in file $result_file does not match that in $header_src\n"
                . "$header_src:" . join(",", @header_fields). "\n"
                . "$result_file:$header_line"
                );
        }
    } else {
        $self->_header_fields(\@header_fields);
        $self->_header_source_file($result_file);
    }

    my @lines;
    while(my $line = $fh->getline) {
        chomp $line;
        push @lines, $line if $line;
    }
    return @lines;
}

sub execute {
    my $self = shift;

    my $json = new JSON;
    my @params = map {$json->decode($_)} $self->job_params;
    my %phenotypes = map {$_->{phenotype} => 1} @params;
    my %genes = map {$_->{gene} => 1} @params;

    my @lines;
    my @null_pairs;
    for my $phenotype (sort keys %phenotypes) {
        for my $gene (sort keys %genes) {
            #for each gene, find files, strip header and then aggregate result lines
            #if there is a null file then go ahead and enter a NULL line
            my $file_stem = $self->results_directory . "/${phenotype}_${gene}";
            my $result_file = "$file_stem.burden.csv";
            my $null_file = "$file_stem.null";
            my $error_file = "$file_stem.error";
            if(-e $result_file) {
                push @lines, $self->_parse_result_file($result_file);
            }
            elsif(-e $null_file) {
                push @null_pairs, [$phenotype, $gene];
            }
            elsif(-e $error_file) {
                $self->error_message("Error calculating burden test for $phenotype and gene $gene")
            }
            else {
                confess $self->error_message(
                    "None of the expected output files for $phenotype and gene $gene produced");
            }
        }
    }

    my $ofh = Genome::Sys->open_file_for_overwriting($self->output_file);
    if (@lines) {
        my $header_fields = $self->_header_fields;
        my $num_fields = scalar @$header_fields;
        my $header = join(",", @$header_fields);
        $ofh->print("$header\n");
        $ofh->print(join("\n",
            @lines, # non-NULL pheno/gene pairs
            (map { join(",", @$_, ("NA") x ($num_fields-2)) } @null_pairs) # NULL pheno/gene pairs
            ) . "\n"
        );
    }

    return 1;
}

1;
