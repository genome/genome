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
            default_value => 'id,hostname,is_current,is_attached,attached_via,desc,ftp_detail,nfs_detail,base_dir'
        },
    ],
    doc => 'list known Genome Modeling Systems'
};

1;

