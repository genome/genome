package Genome::Config::AnalysisProject::SubjectMapping::Command::View;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text qw(justify);

class Genome::Config::AnalysisProject::SubjectMapping::Command::View {
    is => [
        'Genome::Command::Viewer',
        'Genome::Command::WithColor',
    ],
    has_input => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'the analysis project that owns the subject mappings to view',
            shell_args_position => 1,
        },
    ],
    
    has_optional_input => [
        COLUMN_WIDTH => {
            is => 'Number',
            default_value => 40,
        },

    ],
};

sub write_report {
    my ($self, $width, $handle) = @_;
    $self->_write_heading($handle);
    $self->_write_subject_mappings($handle);

    1;
}

sub _write_heading {
    my ($self, $handle) = @_;
    $self->_write_analysis_project_summary($handle, $self->analysis_project);
}

sub _write_analysis_project_summary {
    my ($self, $handle, $analysis_project) = @_;
    $self->_write_analysis_project_heading($handle, $analysis_project);
    $self->_write_analysis_project_body($handle, $analysis_project);
}

sub _write_analysis_project_heading {
    my ($self, $handle, $analysis_project) = @_;
    print $handle 'Analysis Project: ', $analysis_project->id, "\n";
}

sub _write_analysis_project_body {
    my ($self, $handle, $analysis_project) = @_;
    $self->_write_analysis_project($handle, $analysis_project);
    print $handle "\n";
}

sub _write_analysis_project {
    my ($self, $handle, $analysis_project) = @_;
    for my $line ($self->_get_analysis_project_lines($analysis_project)) {
        $self->_write_pairs_line($handle, @$line);
    }
}

sub _write_subject_mappings {
    my ($self, $handle) = @_;

    for my $subject_mapping ($self->analysis_project->subject_mappings){
        $self->_write_subject_mapping_heading($handle, $subject_mapping);
        $self->_write_subject_mapping_body($handle, $subject_mapping);
    }
}

sub _write_subject_mapping_heading {
    my ($self, $handle, $subject_mapping) = @_;
    print $handle 'Subject Mapping: ', $subject_mapping->id, "\n";

}

sub _write_subject_mapping_body {
    my ($self, $handle, $subject_mapping) = @_;
    $self->_write_subjects($handle, $subject_mapping);
}

sub _write_subjects {
    my ($self, $handle, $subject_mapping) = @_;

    for my $subject ($subject_mapping->subject_bridges){
        $self->_write_subject_heading($handle, $subject);
        $self->_write_subject_body($handle, $subject);
        print $handle "\n";
    }
}

sub _write_subject_heading {
    my ($self, $handle, $subject) = @_;
    print $handle '    ';
    print $handle 'Subject: ', $subject->id, "\n";
}

sub _write_subject_body {
    my ($self, $handle, $subject) = @_;

    for my $line ($self->_get_subject_lines($subject)) {
        print $handle '        ';
        $self->_write_pairs_line($handle, @$line);
    }
}

sub _get_subject_lines {
    my ($self, $subject) = @_;
    return (
        ['ID', $subject->id, ],
        ['Subject', $subject->subject->id ],
        ['Label', $subject->label],
        ['Created By', $subject->created_by], 
        ['Created At', $subject->created_at], 
        ['Updated At', $subject->updated_at], 
        ['Tags', join(',', map{$_->name} $subject->subject_mapping->tags)],
    );
}

sub _get_analysis_project_lines {
    my ($self, $analysis_project) = @_;
    return (
        ['Name', $analysis_project->name],
        ['Status', $analysis_project->status],
        ['Run as', $analysis_project->run_as],
        ['Created At', $analysis_project->created_at],
        ['Updated', $analysis_project->updated_at],
        ['Created By', $analysis_project->created_by],
    );
}

sub _write_pairs_line {
    my ($self, $handle, $l_label, $l_value, $r_label, $r_value) = @_;

    if ($r_label and $r_value) {
        print $handle justify($self->_color_pair($l_label, $l_value), 'left',
            $self->COLUMN_WIDTH), " ", $self->_color_pair($r_label, $r_value), "\n";

    } else {
        print $handle $self->_color_pair($l_label, $l_value), "\n";
    }
}

1;
