package Genome::Model::Tools::Vcf::FinalPassOnly;

use Genome;
use strict;
use warnings;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;

class Genome::Model::Tools::Vcf::FinalPassOnly {
    is => 'Command::V2',
    has => [
        output_file => {
            is  => 'Text',
            doc => "output vcf file, it can be either .vcf or .vcf.gz",
        },
        input_file => {
            is  => 'Text',
            doc => "input vcf file, it can be either .vcf or .vcf.gz",
        },
        sample_name => {
            is  => 'Text',
            doc => "Sample name, usually it is tumor sample name",
        },
        tabix_index => {
            is  => 'Boolean',
            doc => 'create tabix index for output vcf file',
            is_optional => 1,
            default     => 1,
        },
    ],
};


sub help_detail {
    <<'HELP';
   Make a new vcf file that only contains final passed variants

HELP
}


sub execute {
    my $self = shift;
    my $out_file = $self->output_file;

    my $in_vcf  = Genome::File::Vcf::Reader->new($self->input_file);
    my $header  = $in_vcf->header;
    my $out_vcf = Genome::File::Vcf::Writer->new($out_file, $header);

    my $index = $header->index_for_sample_name($self->sample_name);

    while (my $entry = $in_vcf->next) {
        next if $entry->is_filtered;
        my $filters = $entry->filters;
        if ($entry->sample_field($index, 'FT') eq 'PASS' or $filters ~~ ['PASS']) {
            $out_vcf->write($entry);
        }
    }
    $out_vcf->close;

    if ($self->tabix_index) {
        if (Genome::Sys->file_is_gzipped($out_file)) {
            my $cmd = Genome::Model::Tools::Tabix::Index->create(
                input_file => $out_file,
                preset => 'vcf',
            );
            unless ($cmd->execute) {
                $self->fatal_message("Could not tabix index $out_file");
            }
        }
        else {
            $self->warning_message("output file $out_file is not a compressed file that tabix requires.");
        }
    }
    return 1;
}

