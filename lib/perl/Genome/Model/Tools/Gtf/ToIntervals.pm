package Genome::Model::Tools::Gtf::ToIntervals;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::ToIntervals {
    is => ['Command'],
    has_input => [
        gtf_file => {
            is => 'Text',
            doc => 'The path to a GFF/GTF format file with gene_id and transcript_id attributes.',
        },
        seqdict_file => {
            is => 'Text',
            doc => 'The sequence dictionary for a reference genome of which the annotation is based on.',
        },
        interval_file => {
            is => 'Text',
            doc => 'The output Interval format file.',
            is_output => 1
        },
    ],
};

sub help_brief {
    "Convert a GTF file into a Interval file.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt ref-cov gtf-to-interval...
EOS

};

sub help_detail {
    'This command will take a GTF format file and convert it to a Interval format file that contains a name(field 5) value that matches RefCov expected values for generating merged metrics.  The format of the name is:
$GENE_ID:$TRANSCRIPT_ID:$TYPE:$ORDINAL:$DIRECTION
'
}

sub execute {
    my $self = shift;

    my $gff_reader = Genome::Utility::IO::GffReader->create(
        input => $self->gtf_file,
    );
    my @interval_headers = qw/chr start end strand name/;
    my $tmp_file = Genome::Sys->create_temp_file_path();
    my $interval_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $tmp_file,
        headers => \@interval_headers,
        separator => "\t",
        print_headers => 0,
        ignore_extra_columns => 1,
    );
    while (my $data = $gff_reader->next_with_attributes_hash_ref) {
        my $attributes = delete($data->{attributes_hash_ref});
        my $ordinal = 'na';
        if ($attributes->{exon_number}) {
            $ordinal = $attributes->{exon_number};
            # This seems like a bug with the parser...
            $ordinal =~ s/\"//g;
        }
        $data->{name} = $attributes->{gene_id} .':'. $attributes->{transcript_id} .':'. $data->{type} .':'. $ordinal .':'. $data->{strand};
        $interval_writer->write_one($data);
    }
    $interval_writer->output->close;
    
    my $cmd = 'cat '. $self->seqdict_file .' '. $tmp_file .' > '. $self->interval_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->seqdict_file,$tmp_file],
        output_files => [$self->interval_file],
    );
    return 1;
}


1;
