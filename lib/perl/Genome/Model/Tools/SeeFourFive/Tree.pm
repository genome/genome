package Genome::Model::Tools::SeeFourFive::Tree;

use strict;
use warnings;
use Genome;            
use IO::File;
use GraphViz;

class Genome::Model::Tools::SeeFourFive::Tree {
    is => 'Command',
    has => [ 
    c45_file => {},
    lines           => {},
    subtree_hashref_of_arrayrefs => {is_optional=>1},
    debug_mode => {default=>0},
    probability_mode => {default=>0},
    leaf_count => {default=>0},
    ],
};

sub perl_src {
    my $self = shift;
    my $lines=shift;
    unless($lines) {
        $lines = $self->lines;
    }  
    chomp @$lines;
    unless (@$lines) {
        return "# no lines $self->{n}?";
    }
    my $src;
    my %metric_names;

    #this may be ok, but perhaps a hash of each level of tree to figure out if we printed something
    #at that depth before would allow us to use
    #if/elsif correctly
    my $tree_depth_count=0;
    for my $line (@$lines) {
        #use tr to count how deep we are in the tree and remove the | so our previous
        #regex works
        my $current_depth= $line =~ tr/|/ /;

        #three cases here.  
        #Equal to current depth= if
        #Greater than current depth: if
        #less than current depth: output difference # of clsoe brackets

        #my ($metric, $op, $value, $decision) = ($line =~ /\s+([\w\-]+)\s+(\S+)\s+(\S+)\s+:\s+(\S+)/);
        my ($metric, $op, $value, $decision, $total_trials, $errors) = ($line =~ /\s*([\w\-]+)\s+(\S+)\s+([\S\w]+?)\s*:\s*([^\s]*)\s*\(?([\d.]*)\/?([\d.]*)\)?/);
        unless($metric && $op && (defined $value)) {
            $self->error_message("Parsing failure on line: $line\n");
            print "Metric: $metric\n";
            print "OP: $op\n";
            print "value: $value\n";
            print "total trials: $total_trials\n";
            print "ErrorS: $errors\n"; 
            return;
        }
        if(defined $decision) {
            unless(defined $total_trials && defined $errors) {
                $self->error_message("a decision does not have metrics about its error rate: $line");
                return;
            }
        }
        if ($value =~ m/^\s*\D+\s*$/) {
            $op = 'eq'; 
            $value = "'$value'"; 
        }
        $op = '==' if $op eq '=';

        ##protect us o lord from the evil aliens
        $metric =~ s/\s/SPACE/g;
        $metric =~ s/\+/PLUS/g;
        $metric =~ s/\-/MINUS/g;

        if($current_depth >= $tree_depth_count) {
            $src .= "if(\$$metric $op $value) {\n";
        }elsif($current_depth < $tree_depth_count) {
            while ($tree_depth_count > $current_depth) {
                $src .= "}\n";
                $tree_depth_count--;
            }
            $src .= "if(\$$metric $op $value) {\n";
        }
        #having a [ implies a subtree
        if($decision && ($decision !~ m/\[/)) { 
            unless($self->debug_mode || $self->probability_mode) {
                $src .= "return '$decision';\n}\n";
            }
            else {
                my @things_to_output;

                if($self->probability_mode) {
                    ##Laplace correction: correct classifications + 1 / $total_trials + # of classes(this is hardcoded for now)
                    my $probability = (($total_trials - $errors) + 1) / ($total_trials + 2);
                    push @things_to_output, $probability;
                }
                if($self->debug_mode) {
                    my $leaf_count= $self->leaf_count;
                    $leaf_count++;
                    $self->leaf_count($leaf_count);
                    push @things_to_output, $leaf_count;
                }
                $src .= "return ('$decision', " . join(',', @things_to_output) . ");\n}\n";
            }
            $tree_depth_count--;
        }
        elsif($decision) {
            $decision =~ tr/[]//d;
            my $subtree_lines= $self->subtree_hashref_of_arrayrefs->{$decision};
            my $subtree_src = $self->perl_src($subtree_lines);
            $src .= $subtree_src;
            $src .= "}\n";
        }
        $tree_depth_count=$current_depth;
    }
    while ($tree_depth_count > 0) {
        $src .= "}\n";
        $tree_depth_count--;
    }
    return $src;
}

sub load_trees {
    my $self=shift;
    my $c45_file=$self->c45_file;

    my $fh = IO::File->new($c45_file);
    unless($fh) {
        $self->error_message("Supplied file cannot be found/opened: $c45_file : $!");
        return;
    }

    while(my $line=$fh->getline) {
        if ($line =~ m/^Decision Tree/) {
            #store the main decision tree
            my @decision_tree;
            #first throw away blank line
            $fh->getline;

            while(my $line = $fh->getline){
                last if ($line =~ m/^$/); 
                push @decision_tree, $line;
            }
            $self->lines(\@decision_tree);
        }
        elsif (my ($subtree_name) = ($line =~ m/^Subtree \[(.*)\]/)) {
            my $subtreehash=$self->subtree_hashref_of_arrayrefs;
            my @subtree_lines;
            #throw away blank line!
            $fh->getline;

            while((my $line = $fh->getline)) {
                last if ($line =~ m/^$/); 
                push @subtree_lines, $line;
            }
            $subtreehash->{$subtree_name}=\@subtree_lines;
            $self->subtree_hashref_of_arrayrefs($subtreehash);
        }
    }
    unless($self->lines) {
        return 0;
    }
    return 1;
}


sub generate_callback_for_headers {
    my $self = shift;
    my @headers = @_;
    my $src = $self->tree_evaluation_function(@headers);
    print $src;
    my $callback = eval $src;
    unless ($callback) {
        die $@;
    }
    return $callback;
}

sub tree_evaluation_function {
    my $self = shift;
    my @headers = @_;
    my $src;
    $src .= "sub {\n";
    $src .= "    my \$data = \$_[0];\n";
    $src .= "    my (" 
    . join(",", map { '$' . $_ } @headers ) 
    . ') = @$data{@headers};'
    . "\n";

    $src .= $self->perl_src();
    $src .= "};\n";

    return $src;
}


sub as_graphviz_obj {
    my ($self,$graph,$root,$identifier,$lines) = @_;
    unless($lines) {
        $lines = $self->lines;
    }  
    chomp @$lines;
    unless (@$lines) {
        return "# no lines $self->{n}?";
    }
    my $src;
    my %metric_names;
    unless($graph) {
        $graph = GraphViz->new(ratio => .6);
    }
    $$identifier ||= 0;

    #this may be ok, but perhaps a hash of each level of tree to figure out if we printed something
    #at that depth before would allow us to use
    #if/elsif correctly
    my $tree_depth_count=0;
    my @parents;
    unless($root) {
        $graph->add_node("All Data");
        @parents = ("All Data");
    }
    else {
        @parents = ($root);
    }

    for my $line (@$lines) {
        #use tr to count how deep we are in the tree and remove the | so our previous
        #regex works
        my $current_depth= $line =~ tr/|/ /;

        #three cases here.  
        #Equal to current depth= if
        #Greater than current depth: if
        #less than current depth: output difference # of clsoe brackets

        #my ($metric, $op, $value, $decision) = ($line =~ /\s+([\w\-]+)\s+(\S+)\s+(\S+)\s+:\s+(\S+)/);
        my ($metric, $op, $value, $decision, $total_trials, $errors) = ($line =~ /\s*([\w\-]+)\s+(\S+)\s+([\S\w]+?)\s*:\s*([^\s]*)\s*\(?([\d.]*)\/?([\d.]*)\)?/);
        unless($metric && $op && (defined $value)) {
            $self->error_message("Parsing failure on line: $line\n");
            print "Metric: $metric\n";
            print "OP: $op\n";
            print "value: $value\n";
            print "total trials: $total_trials\n";
            print "ErrorS: $errors\n"; 
            return;
        }
        if(defined $decision) {
            unless(defined $total_trials && defined $errors) {
                $self->error_message("a decision does not have metrics about its error rate: $line");
                return;
            }
        }
        if ($value =~ m/[A|G|C|T]/) {
            $op = 'eq'; 
            $value = "'$value'"; 
        }
        $op = '==' if $op eq '=';

        ##protect us o lord from the evil aliens
        if($current_depth >= $tree_depth_count) {
            #NEW NODE
            $graph->add_node($$identifier, label => "$metric $op $value");
            if(@parents) {
                $graph->add_edge($parents[-1] => $$identifier);
            }
            push @parents,$$identifier;

            ++$$identifier;

        }elsif($current_depth < $tree_depth_count) {
            while ($tree_depth_count > $current_depth) {
                pop @parents;
                $tree_depth_count--;
            }

            $graph->add_node($$identifier, label => "$metric $op $value");
            if(@parents) {
                $graph->add_edge($parents[-1] => $$identifier);
            }
            push @parents, $$identifier;
            ++$$identifier;
        }
        #having a [ implies a subtree
        if($decision && ($decision !~ m/\[/)) { 
#            my $color = ($decision eq 'WT') ? 'red' : 'green'; 
#            $graph->add_node($$identifier, label => "$decision", style => 'filled', fillcolor => $color);
            $graph->add_node($$identifier, label => "$decision",);
            if(@parents) {
                $graph->add_edge($parents[-1] => $$identifier);
                ++$$identifier;
            }
            $tree_depth_count--;
            pop @parents; #b/c parent is leaf
        }
        elsif($decision) {
            $decision =~ tr/[]//d;
            ##$self->error_message( "Unsupported creation of graph with subtrees");
            ##return;
            #add in subtree
            my $subtree_lines= $self->subtree_hashref_of_arrayrefs->{$decision};
            my $temp_g = $self->as_graphviz_obj($graph,$parents[-1],$identifier,$subtree_lines); #should modify by reference
            pop @parents;
        }
        $tree_depth_count=$current_depth;
    }

    return $graph;
}








1;

