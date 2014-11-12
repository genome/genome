package Genome::Command;

use strict;
use warnings;

use Genome;

class Genome::Command {
    is => 'Command::Tree',
};

# This map allows the top-level genome commands to be whatever
# we wish, instead of having to match the directory structure.
my %command_map = (
    'variant-reporting' => 'Genome::VariantReporting::Command',
    'analysis-project' => 'Genome::Config::AnalysisProject::Command',
    'config' => 'Genome::Config::Command',
    'db' => 'Genome::Db::Command',
    'disk' => 'Genome::Disk::Command',
    'druggable-gene' => 'Genome::DruggableGene::Command',
    'feature-list' => 'Genome::FeatureList::Command',
    'individual' => 'Genome::Individual::Command',
    'instrument-data' => 'Genome::InstrumentData::Command',
    'library' => 'Genome::Library::Command',
    'model' => 'Genome::Model::Command',
    'model-group' => 'Genome::ModelGroup::Command',
    'population-group' => 'Genome::PopulationGroup::Command',
    'processing-profile' => 'Genome::ProcessingProfile::Command',
    'project' => 'Genome::Project::Command',
    'project-part' => 'Genome::ProjectPart::Command',
    'report' => 'Genome::Report::Command',
    'sample' => 'Genome::Sample::Command',
    'search' => 'Genome::Command::Search',
    'software-result' => 'Genome::SoftwareResult::Command',
    'subject' => 'Genome::Subject::Command',
    'sys' => 'Genome::Sys::Command',
    'task' => 'Genome::Task::Command',
    'taxon' => 'Genome::Taxon::Command',
    'tools' => 'Genome::Model::Tools',
);

$Genome::Command::SUB_COMMAND_MAPPING = \%command_map;

use Genome::Command::Base;
%Command::V2::ALTERNATE_FROM_CLASS = (
    %Command::V2::ALTERNATE_FROM_CLASS,
    %Genome::Command::Base::ALTERNATE_FROM_CLASS,
);

1;
