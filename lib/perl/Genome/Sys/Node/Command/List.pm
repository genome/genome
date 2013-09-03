#!/usr/bin/env genome-perl
use strict;
use warnings;
use Genome;

package Genome::Sys::Node::Command::List;

class Genome::Sys::Node::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => { 
            is_constant => 1,
            default_value => 'Genome::Sys::Node',
        },
        show => {
            default_value => 'id,hostname,is_current,is_attached,desc,ftp_detail,nfs_detail'
        },
    ],
    doc => 'list known Genome Modeling Systems'
};

1;

