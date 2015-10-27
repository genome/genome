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

    my %position_score;

    my $type      = $self->output_filter;
    my $output_fh = Genome::Sys->open_file_for_writing($self->parsed_file);

    my ($netmhc_results, $epitope_seq) = $self->parse_input;

    $self->print_header($output_fh);

    my $protein_type = 'MT';
    for my $protein_name (sort keys %{$netmhc_results->{$protein_type}}) {
        for my $variant_aa (sort keys %{$netmhc_results->{$protein_type}->{$protein_name}}) {
            %position_score = %{$netmhc_results->{$protein_type}{$protein_name}{$variant_aa}};
            my @positions = sort {$position_score{$a} <=> $position_score{$b}} keys %position_score;
            my $total_positions = scalar @positions;

            if ($type eq 'all') {
                for (my $i = 0; $i < $total_positions; $i++) {

                    if ($epitope_seq->{'MT'}->{$protein_name}->{$variant_aa}->{$positions[$i]} ne
                        $epitope_seq->{'WT'}->{$protein_name}->{$variant_aa}->{$positions[$i]})
                        # Filtering if mutant amino acid present
                    {
                        my $position = $positions[$i];
                        $self->print_output_line($output_fh, $protein_name, $variant_aa, $position, $position_score{$position}, $netmhc_results, $epitope_seq);
                    }
                }
            }
            if ($type eq 'top') {
                my $position = $positions[0];
                $self->print_output_line($output_fh, $protein_name, $variant_aa, $position, $position_score{$position}, $netmhc_results, $epitope_seq);
            }
        }
    }

    return 1;
}

sub netmhc_file_headers {
    return qw(protein_label position peptide score);
}

sub print_header {
    my $self      = shift;
    my $output_fh = shift;

    print $output_fh join("\t",
        "Gene Name",
        "Point Mutation",
        "Sub-peptide Position",
        "MT score",
        "WT score",
        "MT epitope seq",
        "WT epitope seq",
        "Fold change"
    ) . "\n";
}

sub print_output_line {
    my ($self, $output_fh, $protein_name, $variant_aa, $position, $position_score, $netmhc_results, $epitope_seq) = @_;

    print $output_fh join("\t", $protein_name, $variant_aa, $position, $position_score) . "\t";
    print $output_fh $netmhc_results->{'WT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";

    print $output_fh $epitope_seq->{'MT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";
    print $output_fh $epitope_seq->{'WT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";
    my $fold_change =
        $netmhc_results->{'WT'}->{$protein_name}->{$variant_aa}->{$position} / $position_score;
    my $rounded_FC = sprintf("%.3f", $fold_change);
    print $output_fh $rounded_FC . "\n";
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

    my (%netmhc_results, %epitope_seq);
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
        my ($protein_type, $protein_name, $variant_aa) = split(/\./, $protein_identifier);

        $netmhc_results{$protein_type}{$protein_name}{$variant_aa}{$position} = $score;
        $epitope_seq{$protein_type}{$protein_name}{$variant_aa}{$position} = $epitope;
    }

    return \%netmhc_results, \%epitope_seq;
}

1;
