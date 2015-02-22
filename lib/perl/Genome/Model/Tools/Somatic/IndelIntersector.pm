package Genome::Model::Tools::Somatic::IndelIntersector;

use warnings;
use strict;

use Genome;
use IO::File;
use Switch;

class Genome::Model::Tools::Somatic::IndelIntersector{
    is => 'Command',
    has => [
        intersection_output => {
            is => 'Text',
            is_output => 1,
            doc => 'store the merged indels in BED formated file'
        }
    ],
    has_optional => [
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut if the intersection_output is already present. Useful for pipelines.',
        },
        report_output => {
            is => 'String',
            is_optional => 1,
            doc => 'Send report here.',
        },
    ],
    has_many => [
        model => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'enter model name or ID in order to specify the indel output of a particular model'
        },
    ],
    has_many => [
        indel_file => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Use this until model/build support is added.',
        },
    ],
};

sub help_brief {
    "Intersect Many Indel Files",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic indel-intersector --model 1234567 --model 1234567 
EOS
}

sub help_detail {
    return <<EOS 
IndelIntersector combines indel files from disperate inputs to create a single list of indels, with indels that seem to be duplicates, removed.
EOS
}

sub execute {
    my $self = shift;
    my %inputs;
    my %files;
=cut
    if(defined($self->model)){
        my @model_inputs = $self->model;
        my @models;
        for my $m (@model_inputs){
            my @results = Command::V2->resolve_param_value_from_text($m,'Genome::Model');
            unless(scalar(@results)==1){
                $self->error_message("Found multiple results for ".$m);
                die $self->error_message;
            }
            push @models, $results[0];
        }
        unless(scalar(@models)>0){
            $self->error_message("Found no input models.");
            die $self->error_message;
        }
        for my $model (@models){
            my $num;
            switch($num = $self->number_of_columns($model)){
                case 5  {$inputs{$model->genome_model_id}="pindel"}
                case 36 {$inputs{$model->genome_model_id}="samtools"}
                else { die $self->error_message("Found an indel file with ".$num." columns. Not currently handled.")}
            }
        }
    } elsif (defined($self->indel_file)){
=cut
        my @indel_files = $self->indel_file;
        for my $file (@indel_files){
            unless(-s $file){
                $self->error_message("Could not stat input file ".$file);
                die $self->error_message;
            }
            my $num;
            switch($num = $self->number_of_columns($file)){
                case 5  {$files{$file}="pindel"}
                case 6  {$files{$file}="pindel"}
                case 36 {$files{$file}="samtools"}
                else { die $self->error_message("Found an indel file with ".$num." columns. Not currently handled.")}
            }
        }
  #  }
    if(defined($self->skip_if_output_present) && $self->skip_if_output_present == 1)
    {
        if(-e $self->intersection_output){
            $self->error_message("An existing output file was found. Exiting.");
            die $self->error_message;
        }
    }
    my $count=0;
    my $reject_count=0;

    my %results;
    my %master_hash;
    for my $file (sort(keys(%files))){
        print "Processing file ".$file."\n";
        my $fh = IO::File->new($file);
        my $type;
        if($files{$file} eq "pindel"){
            $type = 1;
        }
        while (<$fh>){
            my $line = $_;
            chomp $line;
            my @data = split /\t|\,/,$line;
            my ($chrom, $start, $stop, $ref, $var,$allele1,$allele2,$length1,$length2);
            if(defined($type)){
                ($chrom, $start, $stop, $ref, $var) = @data;
            } else {
                ($chrom, $start, undef, undef, $allele1, $allele2, $length1, $length2) = @data;
                ($ref,$var) = $self->find_alleles($allele1,$allele2,$length1,$length2);
                if($var eq "0"){
                    $stop = $start + length $ref;
                } else {
                    $stop = $start + 1;
                }
            }
            my $merge=undef;
            for my $start_pos (-10..10){
                if(exists($master_hash{$chrom}{$start+$start_pos})){
                    $merge = 1;
                    last;
                }
            }
            unless(defined($chrom) && defined($start) && defined($stop) && defined($ref) && defined($var)){
                $self->error_message("Found an undefined data point, asploding");
                my $thing = join "\t",($chrom,$start,$stop,$ref,$var);
                print $thing."\n";
                die $self->error_message;
            }
            if(defined($merge)){
                $reject_count++;
                if(exists($results{$chrom}{$start})){
                    print "Found a duplicate record at ".$chrom."\t".$start."\n";
                }
                $results{$chrom}{$start} = $stop."\t".$ref."\t".$var;
                next;
            }else{
                $count++;
                if(exists($master_hash{$chrom}{$start})){
                    print "Found a duplicate record at ".$chrom."\t".$start."\n";
                }
                $master_hash{$chrom}{$start}=$stop."\t".$ref."\t".$var;
            }
        }
        $fh->close;
    }

    print "count =  ".$count."\n";
    print "reject=  ".$reject_count."\n";
    my $output = IO::File->new(">".$self->intersection_output);
    for my $chrom_key (sort(keys(%master_hash))){
        for my $pos_key (sort(keys(%{$master_hash{$chrom_key}}))){
            
            my ($stop,$ref,$var) = split /\t/,$master_hash{$chrom_key}{$pos_key};
            my $line = join "\t",($chrom_key,$pos_key,$stop,$ref,$var);
            print $output $line . "\n";
        }
    }
    $output->close;
    if(defined($self->report_output)){
        my $rep_out = IO::File->new(">".$self->report_output);
        for my $chrom_key (sort(keys(%results))){
            for my $pos_key (sort(keys(%{$results{$chrom_key}}))){
                my ($stop,$ref,$var) = split /\t/,$results{$chrom_key}{$pos_key};
                my $line = join "\t",($chrom_key,$pos_key,$stop,$ref,$var);
                print $rep_out $line . "\n";
            }
        }
        $rep_out->close;
    }
    return 1;
}

sub get_indel_file {
    my $self = shift;
    my $model = shift;
    my $indel_file;
    #TODO use model object to get the proper indel file and return it

    
    return $indel_file;
}

sub number_of_columns {
    my $self = shift;
    my $file = shift;
    unless(-s $file){
        $self->error_message("Could not stat ".$file);
        die $self->error_message;
    }
    my $fh = IO::File->new($file);
    unless(defined($fh)){
        $self->error_message("Could nto open file ".$file);
        die $self->error_message;
    }
    my $cols = $fh->getline;
    chomp $cols;
    my @cols = split /\t|\,/,$cols;
    unless(scalar(@cols)>0){
        $self->error_message("Found less than 1 column in the indel input file ".$file);
        die $self->error_message;
    }
    return scalar(@cols)
}

sub find_alleles {
    my $self = shift;
    my $allele1 = shift;
    my $allele2 = shift;
    my $length1 = shift;
    my $length2 = shift;


    if($allele1 eq '*') {
        if($length2 > 0 ) {
            return ("0" , $allele2);
        }
        else {
            return ($allele2, "0");
        }
    }
    else {
         if($length1 > 0) {
            return ("0" , $allele1);
        }
        else {
            return ($allele1, "0");
        }
    }  
} 

1;
