package Genome::Model::Tools::Vcf::CutVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use File::Basename;
use DateTime;
use POSIX;
use Sort::Naturally;

use constant {
    CHROM => 0,
    POS => 1,
    ID => 2,
    REF => 3,
    ALT => 4,
    QUAL => 5,
    FILTER => 6,
    INFO => 7,
    FORMAT => 8,
    SAMPLE => 9,
};

class Genome::Model::Tools::Vcf::CutVcf {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            is_optional => 0,
            doc => "The file to output",
        },
        input_file => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'First file to merge',
        },
        columns_to_remove => {
            is => 'Text',
            is_many => 1,
            doc => 'A list of the column headers of those columns to remove',
        },
    ],
};


sub help_synopsis {
    <<'HELP';
Remove columns from a merged vcf, if an alt is removed, per sample data are adjusted as necessary
HELP
}


sub execute {
    my $self = shift;

    my $ifh = Genome::Sys->open_gzip_file_for_reading($self->input_file);
    my $ofh = Genome::Sys->open_gzip_file_for_writing($self->output_file);

    my $done = 0;
    my $header_line;

    #process the header, being careful to keep the column header line
    while(! $done){
        my $line = $ifh->getline;
        chomp $line;
        if($line =~ m/^##/){
            print $ofh $line."\n";
        } else {
            $header_line = $line;
            $done = 1;
        }
    }
    my @bad_columns = $self->columns_to_remove;
    my %bad_columns;
    my %bad_col_idxs;
    for my $key(@bad_columns){
        $bad_columns{$key}=-1;
    }
    my @column_headers = split /\t/, $header_line;

    #set hash value for bad_headers to the column number of that bad column in the file
    for my $idx (0..$#column_headers){
        if(exists($bad_columns{$column_headers[$idx]})){
            $bad_columns{$column_headers[$idx]} = $idx;
            $bad_col_idxs{$idx} = 1;
        }    
    }
    for my $key(keys(%bad_columns)){
        if($bad_columns{$key} == -1){
            die $self->error_message("Could not find ".$key." in the vcf.\n");
        }
    }
    
    my @new_columns;
    for my $col (@column_headers){
        unless(exists($bad_columns{$col})){
            push @new_columns, $col;
        }
    }
    print $ofh join("\t",@new_columns)."\n";


    #process each line of the input
    LINE: while( my $line = $ifh->getline ){
        chomp $line;
        my @cols = split /\t/, $line;
        my $pos = $cols[POS];
        my $ref = $cols[REF];
        my @alts = split /\,/,$cols[ALT];
        my @bad_data;
        my %bc = %bad_columns;
        for my $bad (keys(%bc)){
            if($cols[$bc{$bad}] eq '.'){
               delete $bc{$bad};
            }
        }
        my $fix_alleles=0;
        my %bad_alleles;

        for my $key ( keys(%bc)){    
            my ($gt) = split /\:/, $cols[$bc{$key}];
            my $a1;
            my $a2;
            if ($gt eq '.') {
                $a1 = $a2 = '.';
            }
            else {
                ($a1,$a2) = split /\//, $gt;
            }
            unless(defined($a2)) {
                $a2 = $a1;
                warn "This sample data doesn't have a properly formed GT: $cols[$bc{$key}] so I'm assuming $a1/$a2 is what you meant\n";
            }
            map { if(($_ ne '.') && ($_ != 0)){ $bad_alleles{$_} = 1; } } ($a1,$a2);
        }
        my $remaining_alts = 0;
        if(keys(%bad_alleles)){
            for my $idx (9..$#cols){
                next if $cols[$idx] eq '.';
                if(exists($bad_col_idxs{$idx})){
                    next;
                }
                my ($gt) = split /\:/, $cols[$idx];
                my $a1;
                my $a2;
                if ($gt eq '.') {
                    $a1 = $a2 = '.';
                }
                else {
                    ($a1,$a2) = split /\//, $gt;
                }
                unless(defined($a2)) {
                    $a2 = $a1;
                    warn "This sample data doesn't have a properly formed GT: $cols[$idx] so I'm assuming $a1/$a2 is what you meant\n";
                }

                for ($a1, $a2){
                    if(($_ ne '.') && ($_ != 0)){
                        $remaining_alts = 1;
                    }
                    if(exists($bad_alleles{$_})){
                        delete $bad_alleles{$_};
                    }
                }
            }
            unless($remaining_alts==1){
                next LINE;
            }
        }
        my @bad_column_counts;
        for my $key (keys(%bad_columns)){
            push(@bad_column_counts,$bad_columns{$key});
        }
        for my $column_number (sort {$b <=> $a} @bad_column_counts) {
            splice(@cols,$column_number,1);
        }

        if(keys(%bad_alleles)){
            my @new_alts;
            for my $idx (0.. $#alts){
                unless(exists($bad_alleles{$idx+1})){
                    push @new_alts, $alts[$idx];
                }
            }
            $cols[4] = join(",",@new_alts);
            my %alt_convert;
            my $count=0;
            for my $idx (0..$#new_alts){
                for my $idxb (0..$#alts){
                    if( $new_alts[$idx] eq $alts[$idxb] ){
                        $alt_convert{$idxb+1} = $idx+1;
                        last;
                    }
                }
            }

            for my $idx (9..$#cols){
                next if $cols[$idx] eq ".";
                $cols[$idx]  = $self->fix_sample_data($cols[$idx], \%alt_convert);
            }
        } 
        print $ofh join("\t",@cols)."\n";
    } 

    $ifh->close;
    $ofh->close;

    return 1;
}

sub fix_sample_data {
    my $self = shift;
    my $sample = shift;
    my $alt_convert = shift;
    my $answer = "*"; 
    my @fields = split /\:/, $sample;
    my @gts = split /\//, $fields[0];
    my @new_gt;
    unless ($fields[0] eq ".") {
        for my $gt (@gts){
            if($gt == 0){
                push @new_gt, $gt;
            } elsif(exists( $alt_convert->{$gt})) {
                push @new_gt, $alt_convert->{$gt};
            } else {
                die $self->error_message("An allele has been removed from ALT which is still referenced!");
            }
        }
        $fields[0] = join("/",@new_gt);
    }
    $answer = join(":", @fields);
    return $answer;
}

1;
