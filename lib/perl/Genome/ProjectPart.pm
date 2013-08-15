package Genome::ProjectPart;

use strict;
use warnings;
use Genome;

class Genome::ProjectPart {
    is => 'Genome::Notable',
    table_name => 'subject.project_part',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        entity_class_name => {
            is => 'Text',
            len => 256,
            column_name => 'PART_CLASS_NAME',
            doc => 'Class name of the object to which this part points',
        },
        entity_class_name_pretty => {
            calculate_from => 'entity_class_name',
            calculate => sub {
                my ($class) = @_;
                my @part = split(/::/,$class);
                if ($part[1] eq 'Model') {
                    my $subclass = $part[-1];
                    $subclass =~ s/([A-Z])/ $1/g;
                    $subclass =~ s/^\s+(.+)/$1/;
                    return $subclass . ' Model';
                } else {
                    my $i = 1;
                    if ($part[$i] =~ /^(Sys)$/) {
                        $i++;
                    }
                    return $part[$i];
                }
            }
        },
        entity_id => {
            is => 'Text',
            len => 64,
            column_name => 'PART_ID',
            doc => 'ID of the object to which this part points',
        },
        entity => {
            is => 'UR::Object',
            id_by => 'entity_id',
            id_class_by => 'entity_class_name',
            doc => 'Actual object this project part represents',
        },
        part => {
            is => 'UR::Object',
            id_by => 'entity_id',
            id_class_by => 'entity_class_name',
            doc => 'Actual object this project part represents',
        },
        project_id => {
            is => 'Text',
            len => 64,
        },
        project => {
            is => 'Genome::Project',
            id_by => 'project_id',
            doc => 'Project of which this is a part',
            constraint_name => 'GPP_GP_FK',
        },
    ],
    has_optional => [
        label => {
            is => 'Text',
            len => 100,
            doc => 'The label for the part.',
        },
        role => {
            is => 'Text',
            len => 100,
            doc => 'The role of the part.',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'Represents a single part of a project',
};

1;
