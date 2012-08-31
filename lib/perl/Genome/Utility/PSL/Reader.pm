# FIXME ebelter
# convert to G:U:IO:Reader to remove G:U:Parser
# used by 
# ./Model/Tools/Blat/Cat.pm
# ./Model/Tools/Blat/PslToLayers.pm
#
package Genome::Utility::PSL::Reader;

use strict;
use warnings;

use Genome;

# This only works for blat output ran with the -noHead option

my @header_fields = qw(matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts qSeq tSeq);

class Genome::Utility::PSL::Reader {
    is => 'Genome::Utility::Parser',
    has => [
            separator => {
                          default_value => "\t",
                      },
            header => {
                       default_value => '0',
                   },
            ],
};

sub header_fields {
    return \@header_fields;
}


1;
