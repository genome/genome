package Genome::Model::SomaticVariation::Command::Compare::Variants;
use strict;
use warnings;
use Genome;

class Genome::Model::SomaticVariation::Command::Compare::Variants {
    is => 'Genome::Model::Command::Compare::Base',
    doc => 'compare variant detection results between two builds',
};

sub help_synopsis {
    my $self = shift;
    my $command_name = $self->command_name;
    return <<EOS
genome model somatic-variation compare variants 139523346 139696451 /tmp/output_dir
EOS
}

sub _compare_files {
    my $self = shift;
    my (
        $old_file,
        $new_file,
        $old_only_results_file,
        $new_only_results_file,
        $common_results_file,
    ) = @_;

    Genome::Model::Tools::Joinx::Intersect->execute(
        input_file_a => $old_file,
        input_file_b => $new_file,
        miss_a_file => $old_only_results_file,
        miss_b_file => $new_only_results_file,
        output_file => $common_results_file,
        use_version => '1.7',
    );
}

1;

