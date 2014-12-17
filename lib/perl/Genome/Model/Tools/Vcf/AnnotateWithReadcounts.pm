package Genome::Model::Tools::Vcf::AnnotateWithReadcounts;

use strict;
use warnings;
use Genome;
use Genome::File::BamReadcount::Reader;
use Genome::File::Vcf::BamReadcountParser;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;

my $RC_TAG = 'BRCT';
my $RC_HEADER = sprintf('<ID=%s,Number=1,Type=String,Description="Bam readcount line, gzipped then base-64 encoded">',$RC_TAG);

class Genome::Model::Tools::Vcf::AnnotateWithReadcounts {
    is => 'Command::V2',
    has_input => [
        vcf_file => {
            is => 'File',
            doc => 'The vcf (gzipped or not) file that is to be annotated.',
        },
        readcount_file_and_sample_name => {
            is => 'File',
            is_many => 1,
            doc => 'The readcount file and the name of the sample-column it should be annotated into.  <filename>:<sample_name>',
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
            my $match_count = 0;
            while (entries_match($readcount_entries{$sample_idx}, $vcf_entry)) {
                if ($match_count > 0) {
                    remove_readcount_from_vcf_entry($vcf_entry, $sample_idx);
                }
                else {
                    add_readcount_to_vcf_entry($readcount_entries{$sample_idx}, $vcf_entry, $sample_idx);
                }
                $match_count++;
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

    my %readcount_readers;
    for my $sample_name (keys %{$self->readcount_filenames}) {
        my $sample_idx = eval{$header->index_for_sample_name($sample_name)};
        my $error = $@;
        if ($error) {
            if ($error =~ /Sample name $sample_name not found in header/) {
                my @sample_names = $header->sample_names;
                push @sample_names, $sample_name;
                $header->sample_names(\@sample_names);
                $sample_idx = $header->index_for_sample_name($sample_name);
            }
            else {
                die $@;
            }
        }
        $readcount_readers{$sample_idx} = Genome::File::BamReadcount::Reader->new(
            $self->readcount_filenames->{$sample_name});
    }
    my $vcf_writer = Genome::File::Vcf::Writer->new($self->output_file, $header);
    return ($vcf_reader, $vcf_writer, \%readcount_readers);
}

sub readcount_filenames {
    my $self = shift;

    my %result;
    for my $value ($self->readcount_file_and_sample_name) {
        my ($filename, $sample_name) = split(/:/, $value);
        $result{$sample_name} = $filename;
    }
    return \%result;
}

sub entries_match {
    my ($readcount_entry, $vcf_entry) = @_;

    if (defined $readcount_entry) {
        my $rc_chrom  = $readcount_entry->chromosome;
        my $rc_pos    = $readcount_entry->position;
        my $vcf_chrom = $vcf_entry->{chrom};
        my $vcf_pos   = $vcf_entry->{position};

        if ($rc_chrom eq $vcf_chrom) {
            if ($rc_pos == $vcf_pos) {
                return 1;
            }
            elsif ($vcf_entry->has_deletion) {
                $rc_pos--;
                return 1 if $rc_pos == $vcf_pos;
            }
        }
    } 
    return 0;
}

sub add_readcount_to_vcf_entry {
    my ($readcount_entry, $vcf_entry, $sample_idx) = @_;

    my $readcount_line = Genome::File::Vcf::BamReadcountParser::encode($readcount_entry->to_string);
    $vcf_entry->set_sample_field($sample_idx, $RC_TAG, $readcount_line);
    return;
}

sub remove_readcount_from_vcf_entry {
    my ($vcf_entry, $sample_idx) = @_;
    $vcf_entry->set_sample_field($sample_idx, $RC_TAG, undef);
    return;
}

1;
