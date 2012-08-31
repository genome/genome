package Genome::Model::MetagenomicShotgun::Command::Status;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::Command::Status {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::MetagenomicShotgun',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'The metagenomic shotgun models.',
        },
    ],
    doc => 'show the status of metagenomic shotgun models',
};

sub execute {
    my $self = shift;

    my @status;
    for my $model ( $self->models ) {
        push @status, $self->_status_for_model($model);
    }

    $self->status_message( join("\n", @status) );

    return 1;
}

sub _status_for_model {
    my ($self, $model) = @_;

    my @sub_model_methods = (qw/ contamination_screen_model metagenomic_nucleotide_model metagenomic_protein_model viral_nucleotide_model viral_protein_model /);

    my $length = (sort { $b <=> $a } map { length } @sub_model_methods )[0] + 2;
    my $spacify = sub{ return $_[0].( ' ' x ($length - length($_[0])) ); };

    my $get_status_for_model = sub{
        my ($model) = @_;
        my $build = $model->latest_build;
        if ( not $build ) {
            return $model->id.' NO_BUILD';
        }
        return join(' ', $model->id, $build->id, $build->status, $build->data_directory);
    };

    my @status = ( 
        $spacify->('MAIN MODEL:').' '.$get_status_for_model->($model),
    );
    for my $sub_model_method ( @sub_model_methods ) {
        my $sub_model = $model->$sub_model_method;
        my $sub_status = $spacify->(ucfirst( join(' ', grep { $_ ne 'model' } split('_', $sub_model_method)) ).':');
        $sub_status .= ' '.$get_status_for_model->($sub_model);
        push @status, $sub_status;
    }

    return join("\n", @status);
}

1;

