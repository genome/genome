package Genome::Model::Tools::EpitopePrediction::ParseNetmhcOutput;

use strict;
use warnings;
use Data::Dumper;
use Genome;
use File::Basename qw/fileparse/;

class Genome::Model::Tools::EpitopePrediction::ParseNetmhcOutput {
    is        => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        netmhc_file => {
            is  => 'Text',
            doc => 'Raw output file from Netmhc',
        },
        output_directory => {
            is => 'Text',
            doc => 'Location of the output',
        },
        output_filter => {
            is => 'Text',
            doc =>
                'Type of epitopes to report in the final output - select \'top\' to report the top epitopes in terms of fold changes,  \'all\' to report all predictions ',
            valid_values => ['top', 'all'],
        },
        key_file => {
            is  => 'Text',
            is_optional => 1,
            doc => 'Key file for lookup of FASTA IDs'
        },
        netmhc_version => {
            is => 'Text',
            doc => 'NetMHC version to use',
            valid_values => ['3.0','3.4'],
            default_value => '3.4',
            is_optional => 1,
        },
    ],
    has_output => [
        parsed_file => {
            is => 'Text',
            doc => 'File to write the parsed output',
            calculate_from => ['output_directory', 'netmhc_file', 'output_filter'],
            calculate => q|
                my ($filename, $directories, $suffix) = File::Basename::fileparse($netmhc_file, '.xls');
                return File::Spec->join($output_directory, "$filename.parsed.$output_filter");
            |,
        },
    ],
};

sub help_brief {
    "FOR NETMHC3.4 : Parses output from NetMHC3.4 for MHC Class I epitope prediction; Uses a special key file that could be generated using gmt epitope-prediction fasta-key --help.
     The parsed TSV file contains predictions for only those epitopes that contain the \"mutant\" SNV ",;
}

sub execute {
    my $self = shift;

    my $type      = $self->output_filter;
    my $output_fh = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->parsed_file,
        separator => "\t",
        headers => [$self->headers],
    );

    my $netmhc_results = $self->parse_input;

    my $protein_type = 'MT';
    PROTEIN:
    for my $protein_name (sort keys %{$netmhc_results->{$protein_type}}) {
        VARIANT_AA:
        for my $variant_aa (sort keys %{$netmhc_results->{$protein_type}->{$protein_name}}) {
            my $mt_position_data = $netmhc_results->{$protein_type}{$protein_name}{$variant_aa};
            my $wt_position_data = $netmhc_results->{'WT'}{$protein_name}{$variant_aa};
            my @positions = sort {$mt_position_data->{$a}->{score} <=> $mt_position_data->{$b}->{score}} keys %{$mt_position_data};

            POSITION:
            for my $position (@positions) {
                my $mt_score            = $mt_position_data->{$position}->{score};
                my $mt_epitope_sequence = $mt_position_data->{$position}->{epitope_sequence};
                my ($wt_score, $wt_epitope_sequence, $fold_change);
                if (defined($wt_position_data) && defined($wt_position_data->{$position})) {
                    $wt_score            = $wt_position_data->{$position}->{score};
                    $wt_epitope_sequence = $wt_position_data->{$position}->{epitope_sequence};
                    $fold_change         = sprintf("%.3f", $wt_score / $mt_score);
                    # Skip if mutant amino acid is not present
                    if ($mt_epitope_sequence eq $wt_epitope_sequence) {
                        next POSITION;
                    }
                }
                else {
                    $wt_score = $wt_epitope_sequence = $fold_change = 'NA';
                }

                my %data = (
                    'Gene Name' => $protein_name,
                    'Mutation' => $variant_aa,
                    'Sub-peptide Position' => $position,
                    'MT score' => $mt_score,
                    'WT score' => $wt_score,
                    'MT epitope seq' => $mt_epitope_sequence,
                    'WT epitope seq' => $wt_epitope_sequence,
                    'Fold change' => $fold_change,
                );
                $output_fh->write_one(\%data);

                if ($type eq 'top') {
                    next VARIANT_AA;
                }
            }
        }
    }

    return 1;
}

sub netmhc_file_headers {
    return qw(protein_label position peptide score);
}

sub headers {
    return (
        'Gene Name',
        'Mutation',
        'Sub-peptide Position',
        'MT score',
        'WT score',
        'MT epitope seq',
        'WT epitope seq',
        'Fold change'
    );
}

sub protein_identifier_for_label {
    my $self = shift;

    my $key_fh = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->key_file,
        separator => "\t",
        headers => [qw(new_name original_name)],
    );

    my %key_hash;
    while (my $keyline = $key_fh->next) {
        #Entry_1	>WT.GSTP1.R187W
        my $new_name = $keyline->{new_name};
        my $original_name = $keyline->{original_name};
        $original_name =~ s/>//g;
        $key_hash{$new_name} = $original_name;
    }

    return \%key_hash;
}

sub parse_input {
    my $self     = shift;

    my $input_fh = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->netmhc_file,
        separator => "\t",
        headers => [$self->netmhc_file_headers],
        ignore_lines_starting_with => 'NetMHC|Protein|(?:^$)',
        allow_extra_columns => 1,
    );

    my $protein_identifier_for_label = $self->protein_identifier_for_label if $self->netmhc_version eq '3.4';

    my %netmhc_results;
    while (my $line = $input_fh->next) {
        my $position      = $line->{position};
        my $score         = $line->{score};
        my $epitope       = $line->{peptide};
        my $protein_label = $line->{protein_label};

        my ($protein_identifier);
        if ($self->netmhc_version eq '3.4') {
            if (exists($protein_identifier_for_label->{$protein_label})) {
                $protein_identifier = $protein_identifier_for_label->{$protein_label};
            }
        }
        elsif ( $self->netmhc_version eq '3.0' )  {
            $protein_identifier = $protein_label;
        }
        my ($protein_type, $protein_name, $variant_aa) = split(/\./, $protein_identifier, 3);

        $netmhc_results{$protein_type}{$protein_name}{$variant_aa}{$position} = {
            score => $score,
            epitope_sequence => $epitope
        };
    }

    return \%netmhc_results;
}

sub best_matches {
    my ($mt_sequence, $wt_position_data) = @_;

    my %match_counts = %{match_counts($mt_sequence, $wt_position_data)};
    my $highest_match_count = (sort {$b <=> $a} values %match_counts)[0];
    my %best_matches;
    while (my ($position, $match_count) = each %match_counts) {
        if ($match_count == $highest_match_count) {
            $best_matches{$position} = $match_count;
        }
    }
    return \%best_matches;
}

sub match_counts {
    my ($mt_sequence, $wt_position_data) = @_;

    my $match_counts;
    for my $position (keys %{$wt_position_data}) {
        $match_counts->{$position} = ($mt_sequence ^ $wt_position_data->{$position}->{epitope_sequence}) =~ tr/\0//;
    }
    return $match_counts;
}

1;
