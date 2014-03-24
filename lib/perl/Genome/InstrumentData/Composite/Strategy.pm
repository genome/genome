package Genome::InstrumentData::Composite::Strategy;

use strict;
use warnings;

use Genome;
use Parse::RecDescent;
use Storable;

class Genome::InstrumentData::Composite::Strategy {
    is => 'Command::V2',
    has_input => [
        strategy => {
            is => 'Text',
            doc => 'The sequence of alignments and/or filtrations to perform',
        },
    ],
    has_optional_output => [
        tree => {
            is => 'HASH',
            doc => 'The parsed structure generated from the strategy',
        },
    ],
};

sub execute {
    my $self = shift;

    my $str = $self->strategy;
    my $parser = Parse::RecDescent->new($self->grammar())
        or die "Failed to create parser from grammar";

    my $tree = $parser->startrule($str);
    if ( not $tree ) {
        $self->error_message('Failed to parse strategy: '.$self->strategy);
        return;
    }
    $self->tree($tree);

    return $tree;
}

sub grammar {
    return q{
        {
            sub find_root { my $node = shift; if(exists $node->{parent}) { return &find_root($node->{parent}); } else { return $node; } };
        }
        startrule: initial_data operations merge api_version end
            { $return = { data => $item[1], action => $item[2], then => $item[3], api_version => $item[4] }; }
        | initial_data operations api_version end
            { $return = { data => $item[1], action => $item[2], api_version => $item[3] }; }
        | <error>

        end: /^\Z/

        operations: operation "then" operations
            {
                my $new_parents = $item[1];
                my $existing_leaves = $item[3];
                my $new_leaves = [];

                for my $new_parent (@$new_parents) {
                    my $copy_of_leaves = Storable::dclone($existing_leaves);
                    for my $leaf (@$copy_of_leaves) {
                        &find_root($leaf)->{parent} = $new_parent;
                    }
                    push @$new_leaves, @$copy_of_leaves;
                }

                $return = $new_leaves;
            }
        | operation
            { $item[1]; }

        operation: "(" operation_list ")"
            { $item[2]; }
        | alignment
            { [ $item[1] ]; }
        | filtration
            { [$item[1] ]; }

        operation_list: operations "and" operation_list
            { [@{ $item[1] }, @{ $item[3] }];  }
        | operations
            { $item[1]; }
        alignment: "aligned" "to" reference "and" annotation "using" aligner
            { $return = $item[7]; $item[7]->{reference} = $item[3]; $item[7]->{annotation} = $item[5]; $item[7]->{type} = 'align'; }
        | "aligned" "to" reference "using" aligner
            { $return = $item[5]; $item[5]->{reference} = $item[3]; $item[5]->{type} = 'align'; }

        filtration: "filtered" "using" filter
            { $return = $item[3]; $item[3]->{type} = 'filter'; }

        merge: "then" "merged" "using" merger deduplicate(?) refine(?)
            {
                $return = $item[4]; $item[4]->{type} = 'merge';
                my $leaf = $item[4];
                if(@{$item[5]}) {
                    my $next = $item[5][0];
                    $leaf->{then} = $next; $next->{type} = 'deduplicate';
                    $leaf = $next;
                }

                if(@{$item[6]}) {
                    my $next = $item[6][0];
                    $leaf->{then} = $next; $next->{type} = 'refine';
                    $leaf = $item[6];
                }
            }

        deduplicate: "then" "deduplicated" "using" deduplicator
            { $item[4]; }

        refine: "then" "refined" "to" known_sites "using" refiner refine(?)
            {
                $return = $item[6];
                $item[6]->{known_sites} = $item[4];
                my $leaf = $item[6];
                if(@{$item[7]}) {
                    my $next = $item[7][0];
                    $leaf->{then} = $next;
                    $next->{type} = 'refine';
                    $leaf = $next;
                }

            }
        | "then" "refined" "using" refiner refine(?)
            {
                $return = $item[4];
                my $leaf = $item[4];
                if(@{$item[5]}) {
                    my $next = $item[5][0];
                    $leaf->{then} = $next;
                    $next->{type} = 'refine';
                    $leaf = $next;
                }

            }

        api_version: "api" version
            { $item[2]; }

        aligner: program_spec
            { $item[1]; }

        filter: program_spec
            { $item[1]; }

        merger: program_spec
            { $item[1]; }

        deduplicator: program_spec
            { $item[1]; }

        refiner: program_spec
            { $item[1]; }

        reference: name
            { $item[1]; }

        annotation: name
            { $item[1]; }

        known_sites: name
            { $item[1]; }

        initial_data: name
            { $item[1]; }

        program_spec: name version params
            {
                $return = {
                    name => $item[1],
                    version => $item[2],
                    params => $item[3],
                };
            }
        | name version
            {
                $return = {
                    name => $item[1],
                    version => $item[2],
                    params => '',
                };
            }

        word: /([\w\.:-]|\\\\)+/
            { $item[1]; }

        name: word
            { $item[1]; }
        | <error>

        version: word
            { $item[1]; }
        | <error>

        params: {
                my $txt = extract_bracketed($text, '[]');
                $txt =~ s/^\[(.*)\]$/$1/;
                $txt=~ s/\\\\([\[\]])/$1/g;
                $return = $txt;
            }

        list: <matchrule:$arg{rule}> /$arg{sep}/ list[%arg]
            { $return = [ $item[1], @{$item[3]} ] }
        |     <matchrule:$arg{rule}>
            { $return = [ $item[1]] }
    };
}

1;
