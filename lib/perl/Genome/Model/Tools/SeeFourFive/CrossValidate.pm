package Genome::Model::Tools::SeeFourFive::CrossValidate;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Math::Random;
use List::Util qw /shuffle/;

class Genome::Model::Tools::SeeFourFive::CrossValidate {
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
    iterations =>
    {
        type => 'Integer',
        is_optional => 1,
        default => 10,
        doc => 'The number of different trees to create',
    },
    expected_error =>
    {
        is_optional=> 1,
        doc =>"THE ANSWER YOU ARE SEARCHING FOR",
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
    my $data_fh=IO::File->new($self->data_file);
    
    #This may prove imprudent if there is a sufficiently large data file to read
    my @data_array = $data_fh->getlines;
    @data_array = shuffle(@data_array);    
    # num lines / num iterations to make the required number of trees for testing!!!!
    my $divisor = scalar(@data_array) / $self->iterations;
    
    #create the datasets
    my $tree_instance = 0;
    my $total = 0;
    my $fp = 0;
    my $fn = 0;
    while($tree_instance < $self->iterations) { 
        my @choppable_data=@data_array;
        my $start = $tree_instance * $divisor;
        my @test_set = splice(@choppable_data, $start,$divisor);
     
        #open a handle
        my $output_fh = IO::File->new(">/tmp/itr$tree_instance.test");
        unless(defined($output_fh)) {
            $self->error_message("Couldn't open /tmp/itr$tree_instance.data for writing");
            return;
        }
        #write the file
        for my $line (@test_set) {
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
        my $data_output_fh = IO::File->new(">/tmp/itr$tree_instance.data");
        for my $line (@choppable_data) {
            print $data_output_fh $line;
        }
        $data_output_fh->close;
        #train a tree and save the output
        `c4.5 -u -f /tmp/itr$tree_instance > /tmp/itr$tree_instance.txt`;
        my @numbers = $self->parse_out_error_numbers("/tmp/itr$tree_instance.txt");
        $total += shift @numbers;
        $fp += shift @numbers;
        $fn += shift @numbers;
         
        #generate an fref and add it to the array of tree frefs
        $tree_instance++;
    }
    $total /= $self->iterations;
    $fp /= $self->iterations;
    $fn /= $self->iterations;
    printf "Total Error: %0.2f%% False Positive: %0.2f%% False Negative: %0.2f%%\n",$total,$fp*100, $fn*100; 
    return ($total, $fp, $fn);
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


sub parse_out_error_numbers {
    my $self=shift;
    my $tree_file=shift;
    my $fh = IO::File->new("$tree_file");
    unless($fh) {
        $self->error_message("Unable to open file $tree_file!!!!!!");
        return
    }
    while(my $line = $fh->getline) {
        unless($line =~ m/Evaluation on test data/) {
            next;
        }
        #if we're here we just got this line^^^^

        for( my $i=0; $i<5; $i++){ 
            $fh->getline;
        }

        $line= $fh->getline;
        my ($total_error)= ($line  =~ m/\(\s*([\d.]+)\%\)\s+\(/);

        for( my $i=0; $i<4; $i++){ 
            $fh->getline;
        }
        $line= $fh->getline;

        my ($right, $wrong) = ($line =~ m/^\s+(\d+)\s+(\d+)/);
        my $fp = $wrong / ($wrong + $right);
        $line= $fh->getline;
        ($wrong, $right) = ($line =~ m/^\s+(\d+)\s+(\d+)/);
        unless(defined $wrong && defined $right) {
            $self->warning_message("A file had only one category(fp or fn will be incorrect): $tree_file")
        }
        my $fn = $wrong / ($wrong + $right);
        return ($total_error, $fp, $fn); 
    }
}
