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
        },
        readcount_files => {
            is => 'File',
            is_many => 1,
        },
        sample_names => {
            is => 'Text',
            is_many => 1,
        },
        output_file => {
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my ($vcf_reader, $vcf_writer, $readcount_readers) = $self->get_file_objects();

    my %readcount_entries;
    for my $sample_name ($self->sample_names) {
        $readcount_entries{$sample_name} = $readcount_readers->{$sample_name}->next;
    }

    while (my $vcf_entry = $vcf_reader->next) {
        $vcf_entry->add_format_field($RC_TAG);
        for my $sample_name ($self->sample_names) {
            my $readcount_entry = $readcount_entries{$sample_name};
            if (entries_match($readcount_entry, $vcf_entry)) {
                add_readcount_to_vcf_entry($readcount_entry, $vcf_entry, $sample_name);
                $readcount_entries{$sample_name} = $readcount_readers->{$sample_name}->next;
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
    my $counter = 0;
    for my $sample_name ($self->sample_names) {
        $readcount_readers{$sample_name} = Genome::File::BamReadcount::Reader->new(
            $self->readcount_file_for_index($counter));
        $counter++;
    }
    return ($vcf_reader, $vcf_writer, \%readcount_readers);
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
    my ($readcount_entry, $vcf_entry, $sample_name) = @_;

    my $sample_index = $vcf_entry->{header}->index_for_sample_name($sample_name);
    my $readcount_line = process_readcount_line($readcount_entry->to_string);
    $vcf_entry->set_sample_field($sample_index, $RC_TAG, $readcount_line);
    return;
}

sub process_readcount_line {
    my $line = shift;
    my $zipped;
    gzip \$line => \$zipped, Minimal => 1;
    return encode_base64($zipped, '');
}

sub readcount_file_for_index {
    my $self = shift;
    my $index = shift;
    my @files = $self->readcount_files;
    return $files[$index];
}

1;

