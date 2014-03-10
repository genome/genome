package Genome::Model::Tools::Gtf::ToJunc;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::ToJunc {
    is => ['Command'],
    has_input => [
        gtf_file => {
            is => 'Text',
            doc => 'The path to a GFF/GTF format file with gene_id and transcript_id attributes.',
        },
        junc_file => {
            is => 'Text',
            doc => 'The output \'junc\' format file.',
        },
        gene_label => {
            is => 'Text',
            valid_values => ['id','name'],
        },
    ],
};

sub help_brief {
    "Convert a GTF file into a \'junc\' format file.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt gtf to-junc...
EOS

};

sub help_detail {
    'This command will take a GTF format file and convert it to a \'junc\' format file.
'
}

sub execute {
    my $self = shift;

    my $gff_reader = Genome::Utility::IO::GffReader->create(
        input => $self->gtf_file,
    );

    # Score here is the transcript_count
    my @junc_headers = qw/name chr start end strand score gene_label transcript_id/;
    my $junc_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->junc_file,
        headers => \@junc_headers,
        separator => "\t",
        print_headers => 0,
    );
    
    my $gene_label_key = 'gene_'. $self->gene_label;
    my %transcript_id_to_gene_label;
    $self->debug_message('Mapping transcript_id to '. $gene_label_key .' from GTF file \''. $self->gtf_file .'\'');
    while (my $data = $gff_reader->next_with_attributes_hash_ref) {
        my $attributes = delete($data->{attributes_hash_ref});
        my $transcript_id = $attributes->{transcript_id};
        my $gene_label = $attributes->{$gene_label_key};
        if ( $transcript_id_to_gene_label{$transcript_id} ) {
            if ($gene_label eq $transcript_id_to_gene_label{$transcript_id}) {
                next;
            } else {
                die('Multiple gene labels '. $gene_label .' and '. $transcript_id_to_gene_label{$transcript_id} .' found for transcript: '. $transcript_id);
            }
        } else {
            $transcript_id_to_gene_label{$transcript_id} = $gene_label;
        }
    }
    $gff_reader->input->close;
    
    my $tmp_bed12_file = Genome::Sys->create_temp_file_path('gtf_to_bed12.bed');
    $self->debug_message('Converting GTF format file \''. $self->gtf_file .'\' to BED12 format file \''. $tmp_bed12_file .'\'');
    my $gtf_to_bed12_cmd = Genome::Model::Tools::Gtf::ToBed12->execute(
        bed12_file => $tmp_bed12_file,
        gtf_file => $self->gtf_file,
    );
    my $bed12_reader = Genome::Utility::IO::BedReader->create(
        input => $tmp_bed12_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    my %junc;
    $self->debug_message('Generating junctions from BED12 format file \''. $tmp_bed12_file .'\'');
    while (my $bed12_data = $bed12_reader->next) {
        unless (
            ($bed12_data->{'chrom'} =~ /\w+/) &&
                ($bed12_data->{'chromStart'} =~ /\d+/) &&
                    ($bed12_data->{'chromEnd'} =~ /\d+/)
                ) {
            next;
        }

        my $blocks = $bed12_data->{'blockCount'};
        if ($blocks == 1) { next; }

        my $start = $bed12_data->{'chromStart'};
        my @b_sizes = split(/,/, $bed12_data->{'blockSizes'});
        my @b_offsets = split(/,/, $bed12_data->{'blockStarts'});
        for (my $b = 1; $b < $blocks; $b++) {
            my $chr = $bed12_data->{chrom};
            my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
            my $right = $start + $b_offsets[$b] + 1;
            my $strand = $bed12_data->{strand};
            my $j_name = "$chr:$left-$right";
            my $j_id = "$chr:$left-$right($strand)";
            my $transcript_id = $bed12_data->{name};
            if ($junc{$j_id}) {
                push @{$junc{$j_id}->{transcript_id}}, $transcript_id;
            } else {
                my %junc_data = (
                    'chr' => $chr,
                    'start' => $left,
                    'end' => $right,
                    'name' => $j_id,
                    'strand' => $strand,
                    'transcript_id' => [$transcript_id],
                );
                $junc{$j_id} = \%junc_data;
            }
        }
    }
    for my $junction_id (sort {$a cmp $b} keys %junc) {
        my $transcript_ids = delete($junc{$junction_id}->{transcript_id});
        my $transcript_count = scalar(@$transcript_ids);
        for my $transcript_id (@$transcript_ids) {
            my $gene_label = $transcript_id_to_gene_label{$transcript_id};
            unless ($gene_label) {
                die('Missing '. $gene_label_key .' for transcript: '. $transcript_id);
            }
            my $junc_data = $junc{$junction_id};
            $junc_data->{score} = $transcript_count;
            $junc_data->{gene_label} = $gene_label;
            $junc_data->{transcript_id} = $transcript_id;
            $junc_writer->write_one($junc_data);
        }
    }
    
    return 1;
}


1;
