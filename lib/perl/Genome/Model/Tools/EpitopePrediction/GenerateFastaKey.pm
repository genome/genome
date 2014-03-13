package Genome::Model::Tools::EpitopePrediction::GenerateFastaKey;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::EpitopePrediction::GenerateFastaKey {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    doc => "Generates a key file to lookup original protein names in the output file of NetMHC 3.4 from the original 21-mer FASTA file for wildtype(WT) and mutant(MT) proteins",
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input FASTA file with 21mer sequences for wildtype(WT) and mutant(MT) proteins generated using \'gmt epitope-prediction generat-variant-seq21mer\'',
        },
        output_file => {
            is  => 'Text',
            doc => 'The output Key file for lookup',
        },
    ],
};

sub execute {
    my $self = shift;

    my $input_fh  = Genome::Sys->open_file_for_reading( $self->input_file );
    my $output_fh = Genome::Sys->open_file_for_writing( $self->output_file );

    my $i = 1;
    while ( my $line = $input_fh->getline ) {
        chomp $line;
        if ( $line =~ /^>/ ) {
            my $original_name = $line;
            my $new_name      = "Entry_" . $i;
            print $output_fh join( "\t", $new_name, $original_name ) . "\n";
            $i++;
        }
    }
    return 1;
}

1;
