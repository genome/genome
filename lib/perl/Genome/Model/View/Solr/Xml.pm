package Genome::Model::View::Solr::Xml;

use strict;
use warnings;

use Genome;

use Genome::Utility::List;

class Genome::Model::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my ($model) = @_;
                my $pp_type = $model->type_name();
                my $solr_type;

                if ($model->isa('Genome::Model::ReferenceAlignment') && $model->is_lane_qc) {
                    $solr_type = 'model - lane_qc';
                } elsif ($pp_type eq 'reference alignment') {
                    $solr_type = 'model - alignment';
                } elsif ($pp_type =~ /somatic/i) {
                    $solr_type = 'model - somatic';
                } elsif ($pp_type =~ /rna/i) {
                    $solr_type = 'model - rna';
                } elsif ($pp_type eq 'clin seq') {
                    $solr_type = 'model - ClinSeq';
                } elsif ($pp_type eq 'convergence') {
                    $solr_type = 'model - convergence';
                } elsif ($pp_type =~ /microarray/i) {
                    $solr_type = 'model - microarray';
                } else {
                    $solr_type = 'model - other';
                }
            },
        },
        display_type => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my ($model) = @_;
                my $pp_type = $model->type_name();
                my $solr_type;

                if ($model->isa('Genome::Model::ReferenceAlignment') && $model->is_lane_qc) {
                    $solr_type = 'Lane QC Model';
                } elsif ($pp_type eq 'reference alignment') {
                    $solr_type = 'Alignment Model';
                } elsif ($pp_type =~ /somatic/i) {
                    $pp_type =~ s/\b(\w)/\U$1/g; # capitalize words
                    $solr_type = $pp_type . ' Model';
                } elsif ($pp_type =~ /rna/i) {
                    $solr_type = 'RNA Model';
                } elsif ($pp_type =~ /microarray/i) {
                    $solr_type = 'Microarray Model';
                } else {
                    $solr_type = 'Model';
                }
                return $solr_type;
            },
        },
        last_succeeded_build_id => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $build = $_[0]->last_succeeded_build();
                my $b = $build->id if $build;
            }
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_model_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                return join ('?id=', '/view/genome/model/status.html', $_[0]->genome_model_id());
            },
        },
        display_label1 => {
            is  => 'Text',
            default => 'last build',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $build = $_[0]->last_succeeded_build();
                return 'none' if !$build;
                return join ('?id=', '/view/genome/model/build/status.html',$build->id());
            },
        },
        display_label2 => {
            is  => 'Text',
            default => 'data directory',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $build = $_[0]->last_succeeded_build();
                return 'none' if !$build;
                return Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL}, $build->data_directory());
            },
        },
        display_label3 => {
            is  => 'Text',
            default => 'summary report',
        },
        display_url3 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $build = $_[0]->last_succeeded_build() || return 'none';
                my $data_dir = $build->data_directory() || return 'none';
                my $report_pathname = join('/', $data_dir, 'reports', 'Summary', 'report.html');
                if (! -e $report_pathname) { return 'none'; }
                my $summary = Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL}, $report_pathname);
            },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                 {
                    name => 'creation_date',
                    position => 'timestamp',
                },
                {
                    name => 'build_ids',
                    position => 'content',
                },
                {
                    name => 'sample_names_for_view',
                    position => 'content',
                },
                {
                    name => 'created_by',
                    position => 'content',
                },
                {
                    name => 'run_as',
                    position => 'content',
                },
                {
                    name => 'processing_profile',
                    position => 'content',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name'
                    ]
                },
                {
                    name => 'data_directory',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
                {
                    name => 'subject_name',
                    position => 'content',
                },
                {
                    name => 'individual_common_name',
                    position => 'content',
                }
# Loading instrument data is *so* slow
#                {
#                    name => 'instrument_data',
#                    position => 'content',
#                    perspective => 'default',
#                    toolkit => 'text',
#                    aspects => [
#                        'id',
#                        'run_name',
#                    ]
#                }
            ],
        }
    ]
};

sub display_url {

    my ($self, $i) = @_;
    my $model = $self->subject();
    my $url;

    if ($i == 0) {
        $url = join('?','/view/genome/model/view/status.html',$model->id);
    } elsif ($i == 1) {
        my $build = $model->last_successful_build();
        $url = join('?','/view/genome/model/build/view/status.html',$build->id);
    }

    return $url;
}

#x display_title
#x display_icon_url
#x display_type
#
#display_label1
#display_url1
#
#display_label2
#display_url2
#
#display_label3
#display_url3


