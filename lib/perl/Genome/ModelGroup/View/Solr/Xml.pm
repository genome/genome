package Genome::ModelGroup::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'modelgroup'
        },
        display_type => {
            is  => 'Text',
            default => 'Model Group',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_modelgroup_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { 
                return join ('?id=', '/view/genome/model-group/status.html',$_[0]->id()); 
            },
        },
        display_label1 => {
            is  => 'Text',
            default => 'convergence model',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $cmodel = $_[0]->convergence_model();
                return if !$cmodel;
                return join ('?id=', '/view/genome/model/convergence/status.html',$cmodel->id());
            },
        },
        display_label2 => {
            is  => 'Text',
            default => 'last build',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $cmodel = $_[0]->convergence_model() || return 'none';
                my $build  = $cmodel->last_succeeded_build() || return 'none';
                return join ('?id=', '/view/genome/model/build/convergence/status.html',$build->id());
            },
        },
        display_label3 => {
            is => 'Text',
            default => 'summary report',
        },
        display_url3 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $cmodel = $_[0]->convergence_model() || return 'none';
                my $build  = $cmodel->last_succeeded_build() || return 'none';
                my $data_dir = $build->data_directory() || return 'none';
                my $url = join( '/',
                    $ENV{GENOME_SYS_SERVICES_FILES_URL},
                    $data_dir, 'reports', 'Summary', 'report.html' );
                return $url;
            },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'models',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'subject_name',
                    ],
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
            ]
        }
    ]
};






1;
