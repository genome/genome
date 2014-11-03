package Genome::Model::Tools::Gtf::ToBed12;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::ToBed12 {
    is => ['Command'],
    has_input => [
        gtf_file => {
            is => 'Text',
            doc => 'The path to a GFF/GTF format file with gene_id and transcript_id attributes.',
        },
        bed12_file => {
            is => 'Text',
            doc => 'The BED12 format file',
        },
    ],
};

sub help_brief {
    "Convert a GTF file into a BED12 file.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt gtf to-bed12...
EOS

};

sub help_detail {
    'This command will take a GTF format file and convert it to a BED12 format file.
'
}

sub execute {
    my $self = shift;

    my $tmp_gene_pred = Genome::Sys->create_temp_file_path();
    unless (Genome::Model::Tools::Gtf::ToGenePred->execute(
        input_gtf_file => $self->gtf_file,
        output_file => $tmp_gene_pred,
        extended => 1,
    )->result) {
        $self->error_message('Failed to convert from GTF to genePred!');
    }
    my $gene_pred_reader = Genome::Utility::IO::GenePredReader->create(
        input => $tmp_gene_pred,
    );
    my $bed12_writer = Genome::Utility::IO::BedWriter->create(
        output => $self->bed12_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    while (my $gene_pred_data = $gene_pred_reader->next) {
        my @starts = split(',',$gene_pred_data->{exonStarts});
        my @ends = split(',',$gene_pred_data->{exonEnds});
        my @blockSizes;
        my @blockStarts;
        for (my $i=0; $i < $gene_pred_data->{exonCount}; $i++) {
            push @blockSizes, ($ends[$i] - $starts[$i]);
            # TODO: Are the txStart positions zero-based or 1-based?
            push @blockStarts, ($starts[$i]- $gene_pred_data->{txStart});
        }
        my %bed12_data = (
            chrom => $gene_pred_data->{chrom},
            # TODO: Are the txStart positions zero-based or 1-based?
            chromStart => $gene_pred_data->{txStart},
            chromEnd => $gene_pred_data->{txEnd},
            name => $gene_pred_data->{name},
            score => 1,
            strand => $gene_pred_data->{strand},
            thickStart => $gene_pred_data->{cdsStart},
            thickEnd => $gene_pred_data->{cdsEnd},
            itemRgb => 0,
            blockCount => $gene_pred_data->{exonCount},
            blockSizes => join(',',@blockSizes),
            blockStarts => join(',',@blockStarts),
        );
        $bed12_writer->write_one(\%bed12_data);
    }
    return 1;
}


1;
