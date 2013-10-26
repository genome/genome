package Genome::Model::Tools::Ensembl::TranscriptInfo;

use strict;
use warnings;

use Genome;

my @DEFAULT_HEADERS = qw/
                            gene_id
                            ensembl_gene_id
                            coord_system
                            seq_region_name
                            gene_biotype
                            transcript_id
                            ensembl_transcript_id
                            transcript_biotype
                            gene_name
                        /;


class Genome::Model::Tools::Ensembl::TranscriptInfo {
    is => ['Genome::Model::Tools::Ensembl::Base'],
    has => [
        reference_build_id => {
            is => 'Text',
            doc => 'The build id for the reference genome',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
        transcript_info_file => {
            is => 'Text',
            doc => 'The path to a transcript info file',
        },
    ],
};

sub execute {
    my $self = shift;

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->transcript_info_file,
        separator => "\t",
        headers => \@DEFAULT_HEADERS,
    );
    unless ($writer) {
        die('Failed to load transcript info file for writing: '. $self->transcript_info_file);
    }
    my $reference_build = $self->reference_build;
    my $chromosomes = $reference_build->chromosome_array_ref;
    my $slice_adaptor = $self->slice_adaptor;
    my $data;
    for my $chr (sort @{$chromosomes}) {
        if ($chr =~ /random/) {
            warn('Skipping random chromosome: '. $chr);
            next;
        }
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        unless ($slice) {
            $slice = $slice_adaptor->fetch_by_region( 'supercontig', $chr );
            unless ($slice) {
                $self->error_message('Slice not found for: '. $chr);
                die($self->error_message);
            }
        }
        my @genes = @{ $slice->get_all_Genes() };
        for my $gene (@genes){
            #Get info on the gene object
            $data->{'gene_id'}++;
            $data->{'ensembl_gene_id'} = $gene->stable_id();
            $data->{'gene_biotype'} = $gene->biotype();
            $data->{'coord_system'}  = $slice->coord_system()->name();
            $data->{'seq_region_name'} = $slice->seq_region_name();
            $data->{'gene_name'} = $gene->external_name();
            unless ($gene->is_known()){
                $data->{'gene_name'} = "Unknown";
            }

            #Get the transcripts associated with this gene object
            my %transcripts;
            my @trans_list = @{$gene->get_all_Transcripts()};
            my $trans_count = scalar(@trans_list);
            #Get the biotype for each transcript
            foreach my $transcript (@trans_list){
                $data->{'transcript_id'}++;
                $data->{'ensembl_transcript_id'} = $transcript->stable_id();
                $data->{'transcript_biotype'} = $transcript->biotype();
                $writer->write_one($data);
            }
        }
    }
    return 1;
}


1;
