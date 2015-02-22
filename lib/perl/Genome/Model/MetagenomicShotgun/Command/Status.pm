package Genome::Model::MetagenomicShotgun::Command::Status;

use strict;
use warnings;

use Genome;

require Term::ANSIColor;

class Genome::Model::MetagenomicShotgun::Command::Status {
    is => 'Command::V2',
    has => [
        model => {
            is => 'Genome::Model::MetagenomicShotgun',
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'The metagenomic shotgun model.',
        },
    ],
    doc => 'Status the details of a metagenomic shotgun model',
};


sub help_brief {
    return 'Detailed status and info about a metagenomic shotgun model';
}

sub help_detail {
    return 'This command will list the sub models of metagenomic shotgun model. The listing includes status, pprocessing profile and reference information.';
}

sub execute {
    my $self = shift;

    my $show_main_model = $self->_show_main_model;
    return if not $show_main_model;

    my $show_sub_models = $self->_show_sub_models;
    return if not $show_sub_models;

    return 1;
}

sub _show_main_model {
    my $self = shift;

    my $model = $self->model;
    my $model_status = $self->_model_status($model, 'main_model');
    $self->status_message($model_status);

    return 1;
}

sub _show_sub_models {
    my $self = shift;

    for my $label ( Genome::Model::MetagenomicShotgun->sub_model_labels ) {
        my $method = $label.'_model';
        my $model = $self->model->$method;
        my $pp = $model->processing_profile;
        my $ref = $model->reference_sequence_build;
        my $build = $model->latest_build;
        my $model_status = $self->_model_status($model, $label);
        $self->status_message($model_status);
        #$self->status_message(ucfirst(join(' ', split('_', $label))).' model: '.$model->id);
    }

    return 1;
}

my %status_colors = (
    na          => "yellow",
    'new'       => "yellow",
    scheduled   => "yellow",
    running     => "cyan",
    succeeded   => "green",
    failed      => "red",
    unstartable => "red",
    abandoned   => "magenta",
);
sub _model_status {
    my ($self, $model, $label) = @_;

    my @info;
    my $display_name_format = '%-10s %s';

    my $build = $model->latest_build;
    my $build_display = 'NA';
    my $build_status = 'NA';
    if ( $build ) {
        $build_display = sprintf($display_name_format, $build->id, $build->data_directory);
        $build_status = $build->status;
    }
    push @info, [ Status => Term::ANSIColor::colored($build_status, $status_colors{lc($build_status)}) ];
    push @info, [ Build => $build_display ];
    push @info, [ Model => $model->id.' '.$model->name ];
    push @info, [ InstData => join(' ', map { $_->id } $model->instrument_data) || 'NA' ];

    my $processing_profile = $model->processing_profile;
    push @info , [ PP => sprintf($display_name_format, $model->processing_profile->id, $model->processing_profile->name) ];

    my $ref = eval{ $model->reference_sequence_build; };
    if ( $ref ) {
        push @info , [ Ref => sprintf($display_name_format, $ref->id, $ref->name) ];
        push @info , [ Fasta => $ref->full_consensus_path('fa') ];
    }

    if ( eval{ $processing_profile->read_aligner_name; } ) {
        my $aligner = $processing_profile->read_aligner_name;
        $aligner .= ' (v '.$processing_profile->read_aligner_version.')' if $processing_profile->read_aligner_version;
        $aligner .= ' ['.$processing_profile->read_aligner_params.')' if $processing_profile->read_aligner_params;
        push @info , [ Aligner => $aligner ];
    }

    if ( eval{ $processing_profile->read_trimmer_name; } ) {
        my $trimmer = $processing_profile->read_trimmer_name;
        $trimmer .= ' (v '.$processing_profile->read_trimmer_version.')' if $processing_profile->read_trimmer_version;
        $trimmer .= ' ['.$processing_profile->read_trimmer_params.')' if $processing_profile->read_trimmer_params;
        push @info , [ Trimmer => $trimmer ];
    }

    my $status = Term::ANSIColor::colored( join(' ', map { ucfirst } split('_', $label)), 'underline' )."\n";
    my $format = '%10s';
    for my $info ( @info ) {
        $status .= Term::ANSIColor::colored(sprintf($format, ucfirst($info->[0].': ')), 'bold');
        $status .= $info->[1]."\n";
    }
    $status .= "\n";

    return $status;
}

1;

