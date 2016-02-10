package Genome::Model::SingleSampleGenotype::Command::ListOutputFiles;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::ListOutputFiles {
    is => 'Command::V2',
    doc => 'List all haplotype caller output files for one or more builds',
    has_input => {
        builds => {
            is => "Genome::Model::Build::SingleSampleGenotype",
            is_many => 1,
        },
    },
};

sub execute {
    my $self = shift;
    for my $build ($self->builds) {
        for my $r (sort {($a->intervals)[0] cmp ($b->intervals)[0]} $build->haplotype_caller_result) {
            print join("\t", $build->subject->name,
                _interval_string($r),
                _file_path_for_result($r))."\n";
        }
    }
    return 1;
}

sub _interval_string {
    my $result = shift;
    return join(",", $result->intervals);
}

sub _file_path_for_result {
    my $result = shift;
    return File::Spec->join($result->output_dir, $result->_vcf_filename);
}

sub sub_command_category {
    return 'analyst tools';
}

1;
