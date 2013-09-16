package Genome::Model::Tools::Annotate::Adaptor::Vcf;

use strict;
use warnings;
use IO::File;
use Genome;
use Genome::File::Vcf::Reader;

class Genome::Model::Tools::Annotate::Adaptor::Vcf {
    is => 'Genome::Model::Tools::Annotate',
    has => [
    vcf_file => {
        is  => 'String',
        is_input  => 1,
        doc => 'VCF-- Snvs only',
    },
    output_file => {
        is => 'Text',
        is_input => 1,
        is_output => 1,
        is_optional=>1,
        doc => "Store output in the specified file instead of sending it to STDOUT."
    },
    skip_if_output_present => {
        is => 'Boolean',
        is_input => 1,
        is_optional => 1,
        default => 0,
        doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
    },
    ],
};

sub help_brief {
    "Converts gzipped or plaintext SNV-only vcf into CHR POS POS REF VAR. For multiple ALT alleles, one line per ALT will be output",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt annotate adaptor vcf --vcf-file=hooboy --output-file=jazzercise
EOS
}

sub help_detail {
    return <<EOS
    Converts gzipped or plaintext SNV-only vcf into CHR POS POS REF VAR. For multiple ALT alleles, one line per ALT will be output
EOS
}

# For now assume we have a somatic file already made from the normal and tumor bam files... separate tool for this perhaps
sub execute {
    my $self = shift;
    unless (-s $self->vcf_file) {
        $self->error_message("vcf file: " . $self->somatic_file . " does not exist or has no size");
        die;
    }

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }
    my $output_fh;
    if($self->output_file) {
        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    }
    else {
        $output_fh = IO::Handle->new();
        $output_fh->fdopen(fileno(STDOUT),">");
    }
    my $reader = Genome::File::Vcf::Reader->new($self->vcf_file);

    while (my $entry = $reader->next()) {
        my $chr = $entry->{chrom};
        my $pos = $entry->{position};
        my $ref = $entry->{reference_allele};
        my @alts = @{$entry->{alternate_alleles}};
        for my $var (@alts) {
            $output_fh->print("$chr\t$pos\t$pos\t$ref\t$var\n");
        }
    }
    return 1;
}

1;

