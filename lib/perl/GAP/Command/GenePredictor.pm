package GAP::Command::GenePredictor;

use strict;
use warnings;

class GAP::Command::GenePredictor {
    is => ['GAP::Command'],
    has => [
            fasta_file      => { 
                                is  => 'TEXT', 
                                doc => 'single fasta file',
                                is_input => 1
                               },
            bio_seq_feature => { 
                                 is          => 'ARRAY',
                                 is_optional => 1,
                                 doc         => 'array of Bio::Seq::Feature',
                                 is_output => 1
                               },
           ],
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Write a set of fasta files for an assembly";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}
 
1;
