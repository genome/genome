#review/notes
#This file currently needs to exist because workflow can't pass bare_args.  Long-term solution is to refactor `gmt snp intersect` to not take bare_args.

package Genome::Model::Tools::Somatic::FilterCeuYri;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::Somatic::FilterCeuYri {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => "Store variants that pass the filter in this file"
        },
        variant_file => {
            is => 'Text',
            is_input => 1,
            doc => "File of variants to be filtered... in sniper format (chromosome, position, reference, genotype) tab delimited"
        },
        datasource_file => {
            is => 'Text',
            is_optional => 1,
            default => '/gscmnt/834/info/medseq/imported_variants_data/CEU_YRI_all.snps.snpfilter.s',
            doc => "File of source variation to filter out... this defaults to the current CEU and YRI file (combined)"
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
    ],
};

sub help_brief {
    "runs gmt snp intersect on CEU/YRI files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt filter-ceu-yri --output-file filtered.out --variant-file variants.in 
EOS
}

sub help_detail {                           
    return <<EOS 
Filter out CEU and YRI from the variant file... really this is just a wrapper of gmt snp intersect since workflow wont take bare args...
EOS
}

sub execute {
    my $self = shift;

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $command = Genome::Model::Tools::Snp::Intersect->create(
        f1_only_output => $self->output_file,
        file1 => $self->variant_file,
        file2 => $self->datasource_file
    );

    my $return = $command->execute;

    unless ($return == 1) {
        $self->error_message("Intersect command returned $return, expecting 1");
        die;
    }

    return 1;
}


1;
