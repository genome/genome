package Genome::Model::Tools::Vcf::AnnotateWithReadcounts;

use strict;
use warnings;
use Genome;
use Genome::File::BamReadcount::Reader;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;
use IO::Compress::Gzip qw(gzip);
use MIME::Base64 qw(encode_base64);

my $RC_TAG = 'BRCT';
my $RC_HEADER = sprintf('<ID=%s,Number=1,Type=String,Description="Bam readcount line, gzipped then base-64 encoded">',$RC_TAG);

class Genome::Model::Tools::Vcf::AnnotateWithReadcounts {
    is => 'Command::V2',
    has_input => [
        vcf_file => {
            is => 'File',
            doc => 'The vcf (gzipped or not) file that is to be annotated.',
        },
        readcount_file_and_sample_idx => {
            is => 'File',
            is_many => 1,
            doc => 'The readcount file and what sample-column it should be annotated into.  <filename>:<sample_index>',
        },
        output_file => {
            is_output => 1,
            doc => 'The output file, should end in ".vcf" or ".vcf.gz".  If ".vcf.gz" the output will be zipped using bgzip.',
        },
    ],
};

sub execute {
    my $self = shift;

    my ($vcf_reader, $vcf_writer, $readcount_readers) = $self->get_file_objects();

    my %readcount_entries;
    my @sample_idxs = keys %{$readcount_readers};
    for my $sample_idx (@sample_idxs) {
        $readcount_entries{$sample_idx} = $readcount_readers->{$sample_idx}->next;
    }

    while (my $vcf_entry = $vcf_reader->next) {
        $vcf_entry->add_format_field($RC_TAG);
        for my $sample_idx (@sample_idxs) {
            my $readcount_entry = $readcount_entries{$sample_idx};
            if (entries_match($readcount_entry, $vcf_entry)) {
                add_readcount_to_vcf_entry($readcount_entry, $vcf_entry, $sample_idx);
                $readcount_entries{$sample_idx} = $readcount_readers->{$sample_idx}->next;
            }
        }
        $vcf_writer->write($vcf_entry);
    }

    return 1;
}

sub get_file_objects {
    my $self = shift;

    my $vcf_reader = Genome::File::Vcf::Reader->new($self->vcf_file);
    my $header = $vcf_reader->{header};
    $header->add_format_str($RC_HEADER);

    my $vcf_writer = Genome::File::Vcf::Writer->new($self->output_file, $header);
    my %readcount_readers;
    for my $sample_idx (keys %{$self->readcount_filenames}) {
        $readcount_readers{$sample_idx} = Genome::File::BamReadcount::Reader->new(
            $self->readcount_filenames->{$sample_idx});
    }
    return ($vcf_reader, $vcf_writer, \%readcount_readers);
}

sub readcount_filenames {
    my $self = shift;

    my %result;
    for my $value ($self->readcount_file_and_sample_idx) {
        my ($filename, $sample_idx) = split(/:/, $value);
        $result{$sample_idx} = $filename;
    }
    return \%result;
}

sub entries_match {
    my ($readcount_entry, $vcf_entry) = @_;

    if (defined $readcount_entry and
        $readcount_entry->chromosome eq $vcf_entry->{chrom} and
        $readcount_entry->position == $vcf_entry->{position}) {
        return 1;
    } else {
        return 0;
    }
}

sub add_readcount_to_vcf_entry {
    my ($readcount_entry, $vcf_entry, $sample_idx) = @_;

    my $readcount_line = process_readcount_line($readcount_entry->to_string);
    $vcf_entry->set_sample_field($sample_idx, $RC_TAG, $readcount_line);
    return;
}

sub process_readcount_line {
    my $line = shift;
    my $zipped;
    gzip \$line => \$zipped, Minimal => 1;
    return encode_base64($zipped, '');
}

1;
