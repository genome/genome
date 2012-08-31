package Genome::ProjectPart;

use strict;
use warnings;
use Genome;

class Genome::ProjectPart {
    is => 'Genome::Notable',
    id_generator => '-uuid',
    id_by => [
        id => { is => 'Text',
         },
    ],
    has => [
        entity_class_name => { 
            is => 'Text', 
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
            column_name => 'PART_ID',
            doc => 'ID of the object to which this part points',
        },
        entity => {
            is => 'UR::Object',
            id_by => 'entity_id',
            id_class_by => 'entity_class_name',
            doc => 'Actual object this project part represents',
        },
        project_id => { is => 'varchar2', len => 64, column_name => 'PROJECT_ID'},
        project => {
            is => 'Genome::Project',
            id_by => 'project_id',
            doc => 'Project of which this is a part',
        },
    ],
    has_optional => [
        label => { is => 'Text', doc => 'The label for the part.', },
        role => { is => 'Text', doc => 'The role of the part.'  }
    ],
    table_name => 'GENOME_PROJECT_PART',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a single part of a project',
};

1;
