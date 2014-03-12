package Genome::Model::ReferenceAlignment::Command::CheckLaneQc;

class Genome::Model::ReferenceAlignment::Command::CheckLaneQc {
    is => 'Genome::Command::Base',
    doc => 'make sure lane QC exists for the supplied instrument data',
    has => [
        genome_models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            shell_args_position => 1,
            doc => 'models to get or create lane QC models for',
        },
    ],
};


sub execute {
    my $self = shift;
    
    for my $model ($self->genome_models) {
        my @lane_qc_models = $model->get_lane_qc_models;
        if (!@lane_qc_models) {
            $self->print_message('Model (' . $model->__display_name__ . ') is missing lane QC models!');
        }
        else {
            for my $lane_qc_model (@lane_qc_models) {
                my @builds = $lane_qc_model->builds;
                if (@builds) {
                    $status = $lane_qc_model->latest_build->status;
                }
                elsif ($lane_qc_model->build_requested) {
                    $status = 'Queued';
                }
                else {
                    $status = 'Build Requested';
                    $lane_qc_model->build_requested(1);
                }
                $status .= ' ' x (15 - length($status));
                $self->print_message(join "\t", $status, $lane_qc_model->__display_name__);
            }
        }
    }

    return 1;
}


sub print_message {
    my $self = shift;
    my $msg = shift;
    chomp($msg);

    print "$msg\n";

    return 1;
}
