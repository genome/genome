package Genome::Model::Tools::Simulation::HapgenToBed;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
class Genome::Model::Tools::Simulation::HapgenToBed {
    is => 'Command',
    has_optional_input => [
    output_prefix=>
    {
        is=>'Text',
        is_optional=>1,
        is_output=>1,
        default=>"individual",
    },
    ],
    has_input=> [
    input_haps_file=>
    { 
        is=>'Text',
        is_optional=>0,
    },
    input_legend_file=>
    {
        is=>'Text',
        is_optional=>0,
    },
    ],
    
};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}

sub execute {
    $DB::single=1;
    my $self = shift;
    unless(-s $self->input_legend_file) {
        $self->error_message("legend file has no size.");
        return 0;
    }
    unless(-s $self->input_haps_file) {
        $self->error_message("haplotypes file has no size.");
        return 0;
    }
    my @haps = $self->process_legend_file($self->input_legend_file);
    my $haps_fh = Genome::Sys->open_file_for_reading($self->input_haps_file);
    my $line = $haps_fh->getline;
    chomp($line);
    my @fields = split " ", $line;
    my $number_of_individuals = scalar(@fields)/2;
    if($number_of_individuals ==0 ){
        return 1;
    }
    my @fhs = $self->prepare_files_for_writing($number_of_individuals, $self->output_prefix);
    my $count=0;
    do {
        chomp($line);
        my @fields = split /\s+/, $line;
        for(my $i=0; $i < @fields; $i+=2) {
            my $all0 = $haps[$count][$fields[$i]];
            my $all1 = $haps[$count][$fields[$i+1]];
            my $stop = $haps[$count][2];
            my $start = $stop -1;
            unless ($stop) {
             $DB::single=1;
            }
            my $chr = "22"; #FIXME. do something better if we ever include more chroms
            $fhs[$i/2]->print("$chr\t$start\t$stop\t$all0/$all1\n");
        }

            $count++;
    }while($line = $haps_fh->getline);
    return 1;


        

}
    1;


sub process_legend_file { 
    my $self = shift;
    my $file = shift;
    my $fh = Genome::Sys->open_file_for_reading($file);
    $fh->getline; #throw away header
    my @haps;
    my $count =0;
    while(my $line = $fh->getline) {
        chomp($line);
        my ($id, $position, $all0, $all1)  = split /\s+/, $line;
        $haps[$count][0]=$all0;
        $haps[$count][1]=$all1;
        $haps[$count][2]=$position;
        $DB::single=1 unless $position;
        $count++;
    }
    return @haps;
}

sub prepare_files_for_writing {
    my $self = shift;
    my $number_of_individuals = shift;
    my $output_prefix=shift;
    my @filehandles;
    for (my $i=0; $i < $number_of_individuals; $i++) {
        my $fh = IO::File->new($output_prefix . "_$i.bed" , ">");
        $filehandles[$i]=$fh;
    }
    return @filehandles;
}

