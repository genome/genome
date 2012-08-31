package Genome::Model::MutationalSignificance::Command::CreateBamList;

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command::CreateBamList {
    is => ['Command::V2'],
    has_input => [
        somatic_variation_builds => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_many => 1,    
        },
    ],
    has_input_output => [
        bam_list => {
            is => 'Text',},
    ],
};

sub execute {
    my $self = shift;

    my $out_string = "";

    foreach my $build ($self->somatic_variation_builds) {
        $out_string .= $build->tumor_build->model->subject->extraction_label; 
        $out_string .= "\t";
        $out_string .= $build->normal_build->whole_rmdup_bam_file;
        $out_string .= "\t";
        $out_string .= $build->tumor_build->whole_rmdup_bam_file;
        $out_string .= "\n";
    }

    my $fh = IO::File->new($self->bam_list, '>');
    $fh->print($out_string);
    $fh->close;

    $self->status_message('Created BAM list');
    return 1;
}

1;
