package Genome::Disk::Command::Volume::Resize;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Volume::Resize {
    is => 'Command::V2',
    has => [
        volume => {
            is => 'Genome::Disk::Volume',
            shell_args_position => 1,
        },
    ],
    doc => 'Figures out the size of the supplied volume and updates the database accordingly',
};

sub help_detail {
    return 'DFs the given volume and updates the size of the volume in the database';
}

sub execute {
    my $self = shift;
    my $volume = $self->volume;

    my $mount_path = $volume->mount_path;
    my @df_output = qx(df -Pk $mount_path);
    unless (@df_output == 2) {
        $self->error_message('\'df\' output does not match expected pattern, exiting');
        return;
    }
    my $df_total_kb = (split(/\s+/, $df_output[1]))[1];
    my $volume_total_kb = $volume->total_kb;
    my $delta_total_kb = $df_total_kb - $volume_total_kb;

    if ($delta_total_kb == 0) {
        $self->status_message('No resize needed.');
        return 1;
    }

    $self->status_message('');
    $self->warning_message('This is not the best way to resize a volume and will soon be deprecated.');
    $self->warning_message('The best way to resize a volume will be to visit the following URL and edit in LIMS\' web app.');
    $self->warning_message(sprintf('https://imp-lims.gsc.wustl.edu/entity/disk-volume/%s', $volume->id));
    $self->warning_message(sprintf('New total_kb is %s.', $df_total_kb));
    $self->warning_message(sprintf('New unallocated_kb is %s.', $volume->unallocated_kb + $delta_total_kb));
    my $question = "Would you like to resize $mount_path by $delta_total_kb kb? From $volume_total_kb kb to $df_total_kb kb.";
    if($self->_ask_user_question($question) eq 'yes') {
        $self->resize_by_kb($delta_total_kb);
        $self->status_message('Volume has been resized.');
        return 1;
    }
    else {
        $self->status_message('Aborting due to user response.');
        return;
    }

    return 1;
}

sub resize_by_kb {
    my $self = shift;
    my $delta_total_kb = shift;
    my $volume = $self->volume;

    $volume->total_kb($volume->total_kb + $delta_total_kb);
    $volume->unallocated_kb($volume->unallocated_kb + $delta_total_kb);

    return 1;
}

1;

