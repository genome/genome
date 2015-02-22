package Genome::Model::Tools::EpitopePrediction::BindingFilter;
# Based on myParseReportsTopGenesTopAlleleStringent_4.pl by TWYLIE (Feb 2012)

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);

use Genome::Utility::IO::SeparatedValueReader qw();
use Genome::Utility::IO::SeparatedValueWriter qw();

class Genome::Model::Tools::EpitopePrediction::BindingFilter {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        fof_file => {
            is  => 'FilePath',
            doc => 'FOF containing list of parsed epitope files for different allele-length combinations (same sample)',
        },
    ],
    has_output => [
        output_file => {
            is  => 'FilePath',
            doc => 'Output .xls file containing list of filtered epitopes based on binding affinity for each allele-length combination per gene',
        }
    ],
};

sub help_brief {
    "Takes in a FOF with parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinity."
}

sub execute {
    my $self = shift;

    my $fof_fh = Genome::Sys->open_file_for_reading($self->fof_file);
    my $out_fh = Genome::Utility::IO::SeparatedValueWriter->create(
        output    => $self->output_file,
        separator => "\t",
        headers   => [
            'Mode',
            'Sample',
            'Length',
            'Gene Name',
            'Allele',
            'Point Mutation',
            'Sub Peptide Position',
            'MT Score',
            'WT Score',
            'MT Epitope Seq',
            'WT Epitope Seq',
            'Fold Change',
        ],
    );

    my %prediction;
    my $threshold = 500 ;
    while (my $file = $fof_fh->getline) {
        chomp $file;
        my $basename = basename( $file);
        my @f      = split( /\./, $basename );
        my $sample = $f[0];
        $sample =~ s/_netmhc//g;
        my $allele = $f[1];
        my $length = $f[2];

        my $mode = 'filtered';
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            input     => $file,
            separator => "\t",
            headers   => [
                'Gene Name',
                'Point Mutation',
                'Sub Peptide Position',
                'MT Score',
                'WT Score',
                'MT Epitope Seq',
                'WT Epitope Seq',
                'Fold Change',
            ],
        );
        $reader->next; # skip headers
        while (my $gene_name = $reader->next) {
            $gene_name->{Allele} = $allele;
            push( @{ $prediction{$mode}->{$sample}->{$length}->{genes} }, $gene_name );
        }
    }
    close ($fof_fh);

    my %best;
    foreach my $mode (sort keys %prediction) {
        foreach my $sample (sort keys %{ $prediction{$mode} }) {
            foreach my $length (sort keys %{ $prediction{$mode}->{$sample} }) {
                foreach my $gene (sort @{ $prediction{$mode}->{$sample}->{$length}->{genes} }) {
# BEST
                    if ( $best{$sample}->{$gene->{'Gene Name'}}->{SCORE} ) {
                        if ($gene->{'MT Score'} > $best{$sample}->{$gene->{'Gene Name'}}->{SCORE}) {
                            next;
                        }
                        if ($gene->{'MT Score'} < $best{$sample}->{$gene->{'Gene Name'}}->{SCORE}) {
                            $best{$sample}->{$gene->{'Gene Name'}}->{GENES} = [];
                        }
                    }

                    $best{$sample}->{$gene->{'Gene Name'}}->{SCORE} = $gene->{'MT Score'};
                    $gene->{Sample} = $sample;
                    $gene->{Length} = $length;
                    $gene->{Mode}   = $mode;
                    push( @{ $best{$sample}->{$gene->{'Gene Name'}}->{GENES} }, $gene );
                }
            }
        }
    }

# REPORTING
    foreach my $sample (sort keys %best) {
        foreach my $gene (sort keys %{ $best{$sample} }) {
            foreach my $entry (@{ $best{$sample}->{$gene}->{GENES} }) {
                if ($entry->{'MT Score'} < $threshold) {
                    $out_fh->write_one($entry);
                }
            }
        }
    }

    return 1;
}
