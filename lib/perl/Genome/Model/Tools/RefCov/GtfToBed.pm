package Genome::Model::Tools::RefCov::GtfToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::GtfToBed {
    is => ['Command'],
    has_input => [
        gff_file => {
            is => 'Text',
            doc => 'The path to a GFF/GTF format file with gene_id and transcript_id attributes.',
        },
        bed_file => {
            is => 'Text',
            doc => 'The output BED format file.',
            is_output => 1
        },
    ],
};

sub help_brief {
    "Convert a GTF file into a RefCov BED file.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt ref-cov gtf-to-bed...
EOS

};

sub help_detail {
    'This command will take a GTF format file and convert it to a BED format file that contains a name(field 4) value that matches RefCov expected values for generating merged metrics.  The format of the name is:
$GENE_ID:$TRANSCRIPT_ID:$TYPE:$ORDINAL:$DIRECTION
'
}

sub execute {
    my $self = shift;

    my $gff_reader = Genome::Utility::IO::GffReader->create(
        input => $self->gff_file,
    );
    my $bed_writer = Genome::Utility::IO::BedWriter->create(
        output => $self->bed_file,
    );
    while (my $data = $gff_reader->next_with_attributes_hash_ref) {
        my $attributes = delete($data->{attributes_hash_ref});
        my $ordinal = 'na';
        if ($attributes->{exon_number}) {
            $ordinal = $attributes->{exon_number};
            # This seems like a bug with the parser...
            $ordinal =~ s/\"//g;
        }
        # BED uses zero-based start position;
        $data->{start}--;
        $data->{name} = $attributes->{gene_id} .':'. $attributes->{transcript_id} .':'. $data->{type} .':'. $ordinal .':'. $data->{strand};
        $bed_writer->write_one($data);
    }
    return 1;
}


1;
