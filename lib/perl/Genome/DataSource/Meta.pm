package Genome::DataSource::Meta;

# The datasource for metadata describing the tables, columns and foreign
# keys in the target datasource

use strict;
use warnings;

use UR;

UR::Object::Type->define(
    class_name => 'Genome::DataSource::Meta',
    is => ['UR::DataSource::Meta'],
);


1;
