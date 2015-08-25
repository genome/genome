package Genome::CommandLine;

use strict;
use warnings;

use Genome;

class Genome::CommandLine {
    is => 'Command::Tree',
};

sub command_tree_source_classes {
    return ('Genome::Model::Tools','Genome::Command');
}

sub _categorize_sub_commands {
    my $class = shift;

    my @sub_command_classes = $class->sorted_sub_command_classes;
    my %categories;
    my @order;
    for my $sub_command_class (@sub_command_classes) {
        next if $sub_command_class->_is_hidden_in_docs();
        #next unless $sub_command_class =~ /Genome::B/;

        my $category;
        if ($sub_command_class =~ /Genome::Model::Tools/) {
            $category = 'tools';
        }
        else {
            $category = 'system';
        }

        unless (exists $categories{$category}) {
            if ($category) {
                push(@order, $category)
            } else {
                unshift(@order, '');
            }
            $categories{$category} = [];
        }
        push(@{$categories{$category}}, $sub_command_class);
    }

    return (\@order, \%categories);
}

1;
