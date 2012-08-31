package Genome::Model::Tools::SeeFourFive::TreeBag;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Math::Random;

class Genome::Model::Tools::SeeFourFive::TreeBag {
    is => 'Command',
    has => [
    name_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "config file detailing name and range 
        of columns in the two snp files for C4.5",
    },
    data_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "file of data which C4.5 will use to 
        attempt to construct good rules. If not specified
        a default subset will be created at runtime",
    },
    test_file =>
    {  
        type => 'String',
        is_optional => 0,
        doc => "file of data which C4.5 will use to
        test out the rules it decides on. if not specified,
        good/bad concat will be used",
    },
    iterations =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 2,
        doc => 'The number of different trees to create',
    },
    ]
};

#Description
#This module will create a variable number of training sets and trees to be used as a decision tree ensemble

sub execute {
    my $self=shift;
    unless(-f $self->data_file) {
        $self->error_message("Data file is not a file: " . $self->data_file);
        return;
    }
    unless(-f $self->name_file) {
        $self->error_message("Names file is not a file: " . $self->name_file);
        return;
    }
    unless(-f $self->test_file) {
        $self->error_message("Test file is not a file: " . $self->test_file);
        return;
    }
    my $data_fh=IO::File->new($self->data_file);
    #This may prove imprudent if there is a sufficiently large data file to read
    my @data_array = $data_fh->getlines;

    #create the datasets
    my $tree_instance = 0;
    while($tree_instance < $self->iterations) { 
        my @indices = $self->create_random_set_of_indices(scalar(@data_array));
        #open a handle
        my $output_fh = IO::File->new(">/tmp/itr$tree_instance.data");
        unless(defined($output_fh)) {
            $self->error_message("Couldn't open /tmp/itr$tree_instance.data for writing");
            return;
        }
        #write the file
        for my $line (@data_array[@indices]) {
            print $output_fh $line;
        }
        $output_fh->close;

        #create links for names and test
        eval {
            symlink $self->name_file,"/tmp/itr$tree_instance.names";
        };

        if($@) {            
            $self->error_message("Unable to symlink names file: $@");
            return;
        }

        eval {
            symlink $self->test_file,"/tmp/itr$tree_instance.test";
        };

        if($@) {
            $self->error_message("Unable to symlink test file: $@");
            return;
        }
        #train a tree and save the output
        `/gscmnt/sata180/info/medseq/bshore/amll123t98_annotation/R8/R8/Src/c4.5 -f /tmp/itr$tree_instance > /tmp/itr$tree_instance.txt`;

        #generate an fref and add it to the array of tree frefs
        $tree_instance++;
    }


}

1;

sub help_detail {
    "This module will bag multiple trees",
}

#This function returns a random sampling of the array for creation of a training set
sub create_random_set_of_indices {
    my ($self,$data_size) = @_;
    return random_uniform($data_size, 0, $data_size-1); #should be a list of random integers
}
