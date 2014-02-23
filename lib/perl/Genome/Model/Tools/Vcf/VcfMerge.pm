package Genome::Model::Tools::Vcf::VcfMerge;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;

class Genome::Model::Tools::Vcf::VcfMerge {
    is => 'Command',
    has => [
    output_file => {
        is => 'Text',
        is_output => 1,
        is_optional => 0,
        doc => "Output merged VCF",
    },

    vcf_files => {
        is => 'Text',
        is_optional => 0,
        doc => "comma-seperated list of VCF files containing mutations from the same sample",
    },

    source_ids => {
        is => 'Text',
        is_optional => 1,
        doc => "given a comma-separated list of ids used to identify the source of the input VCF files. (i.e. GATK, samtools, varScan), will label the source in the info field",
    },

    merge_filters => {
        is => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => "Keep the filter information from all the inputs, (even though we keep most fields only from first file)",
    },

    keep_all_passing => {
        is => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => "Only active if merge-filters is TRUE. If a position is labeled PASS in any file, mark it as passing in the ouput file (union). Default is to let filtering from any input source override PASS.",
    },

    require_all_passing => {
        is => 'Boolean',
        is_optional => 1, 
        default => 0,
        doc => "require that variants be called and passing in all input files to be labeled PASS (intersect) Default is to let filtering from any input source override PASS.",
    },


    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Merge multiple VCFs - keep the quality scores from files in desc order"
}


sub help_synopsis {
    <<'HELP';
Merge multiple VCFs - keep the FORMAT lines from files in desc order.
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Merge multiple VCFs. For identical calls made by different algorithms, merge them, keeping the FORMAT/scores from the file that is listed first in the vcf_files string.
HELP
}

###############

sub execute {                               # replace with real execution logic.
    my $self = shift;
    my $vcf_files = $self->vcf_files;
    my $source_ids = $self->source_ids;
    my $merge_filters = $self->merge_filters;
    my $keep_all_passing = $self->keep_all_passing;

    my @vcffiles = split(/,/,$vcf_files);
    if (@vcffiles < 1){
        die ("requires multiple VCF files to be input (comma-sep)")
    }

    my @vcfnames;
    if(defined($source_ids)){
        @vcfnames = split(",",$source_ids);
        if (@vcffiles != @vcfnames){
            die ("requires a source id for each input VCF file")
        }
    }
    my %vcf_size;
    for (my  $i=0; $i<scalar(@vcffiles); $i++) {
        my $vcf = $vcffiles[$i];
        chomp(my $lines = `wc -l $vcf | cut -f1 -d' '`);
        $vcf_size{$i}{'lines'}=$lines;
        $vcf_size{$i}{'file'}=$vcf;
        $vcf_size{$i}{'name'}=$vcfnames[$i];
    }
    @vcfnames=();
    @vcffiles=();
    $self->debug_message("Sorting files according to size, loading smallest one into memory first");
    for my $key (sort {$vcf_size{$a}{'lines'} <=> $vcf_size{$b}{'lines'}} keys %vcf_size) {
        $self->debug_message("File: " . $vcf_size{$key}{'name'} . " Size: " . $vcf_size{$key}{'lines'});
        push @vcfnames, $vcf_size{$key}{'name'};
        push @vcffiles, $vcf_size{$key}{'file'};
    }




    my %varHash;
    my %infoHash;
    my %filterHash;
    my %passHash;
    my @header;
    #hash the first file
    my $inFh = IO::File->new( $vcffiles[0] ) || die "can't open file\n";


    while(my $line = $inFh->getline )
    {
        chomp($line);
        if ($line =~ /^\#/){
            if ($line =~ /##INFO=\<ID\=(\w+),/){
                $infoHash{$1} = $line
            }
            if ($line =~ /##FILTER=\<ID\=(\S+),/){
                $filterHash{$1} = $line
            }
            push(@header,$line);
            next;
        } 

        my @col = split("\t",$line);

        my $chr = $col[0];
        #replace X and Y for sorting
        $chr = "23" if $col[0] eq "X";
        $chr = "24" if $col[0] eq "Y";
        $chr = "25" if $col[0] eq "MT";
        my $id = $col[1] . ":" . $col[3] . ":" . $col[4];

        #add source id
        if(defined($source_ids)){
            if($col[7] eq "."){
                $col[7] = "VC=" . $vcfnames[0];
            } else {
                $col[7] = $col[7] . ";VC=" . $vcfnames[0];
            }
        }

        @{$varHash{$chr}{$id}} = @col;

        if ($self->require_all_passing){
            if ($col[6] eq "PASS"){
                $passHash{$chr}{$id} = 1;
            }
        }

    }
    close($inFh);


    my @newInfo;
    my @newFilters;
    if(defined($source_ids)){
        push(@newInfo,"##INFO=<ID=VC,Number=.,Type=String,Description=\"Variant caller\">");
    }
    if ($self->require_all_passing){
        push(@newFilters,"##FILTER=<ID=intersect,Description=\"Removed during intersection\">");
    }
    #add data from subsequent files if data does not exist in first file
    for(my $i=1; $i<@vcffiles; $i++){
        my $vcf_file=$vcffiles[$i];

        chomp(my @new_headers = `grep "^#" $vcf_file`);

        for my $header (@new_headers) {
            if ($header =~ /^#/){
                if ($header =~ /##INFO=\<ID\=(\w+),/){
                    unless (exists($infoHash{$1})){
                        push(@newInfo,$header)
                    }
                }
                if ($header =~ /##FILTER=\<ID\=(\S+),/){
                    unless (exists($filterHash{$1})){
                        push(@newFilters,$header)
                    }
                }                    
                next;
            } 
        }

        my $output_file = $self->output_file;
        if(-e $output_file) {
            unless(unlink $output_file) {
                $self->error_message("merge target already exists, unable to delete, aborting");
                return 0;
            }
        }
        $self->debug_message("Printing header...");
        $self->print_header($output_file, \@header, \@newInfo, \@newFilters);
        

        my $prev_chr;
        $inFh = IO::File->new( $vcf_file ) || die "can't open file\n";
        while(my $line = $inFh->getline ) {
            chomp($line);
            next if $line =~m/^#/;

            my @col = split("\t",$line);
            my $chr = $col[0];
            next if($chr =~ /NT/);

            #replace X and Y for sorting
            $chr = "23" if $col[0] eq "X";
            $chr = "24" if $col[0] eq "Y";
            $chr = "25" if $col[0] eq "MT";

            $prev_chr = $chr unless $prev_chr; ## load prev_chr the first time through to prevent weirdness
            if($prev_chr ne $chr) {
                $self->debug_message("Finished with $prev_chr, printing to file...");
                $self->print_chromosome($prev_chr, $output_file, $varHash{$prev_chr},  \@vcffiles, $passHash{$prev_chr});
                delete $varHash{$prev_chr};
                delete $passHash{$prev_chr};
            };



            my $id =  $col[1] . ":" . $col[3] . ":" . $col[4];

            if(exists($varHash{$chr}{$id})){
                #add source id
                if(defined($source_ids)){
                    @{$varHash{$chr}{$id}}[7] = @{$varHash{$chr}{$id}}[7] . "," . $vcfnames[$i];
                }

                #filter
                if($merge_filters){

                    #union
                    if($keep_all_passing){ 
                        if( ( @{$varHash{$chr}{$id}}[6] eq "PASS" ) || ($col[6] eq "PASS")){
                            @{$varHash{$chr}{$id}}[6] = "PASS";
                        }


                        #overlap
                    } else { 
                        if( @{$varHash{$chr}{$id}}[6] eq "PASS" ){
                            @{$varHash{$chr}{$id}}[6] = $col[6];

                            if ($self->require_all_passing){
                                if ($col[6] eq "PASS"){
                                    $passHash{$chr}{$id} = $passHash{$chr}{$id} + 1;
                                }
                            }

                        } else {
                            unless ( $col[6] eq "PASS" ){
                                @{$varHash{$chr}{$id}}[6] = @{$varHash{$chr}{$id}}[6] . ";" . $col[6];
                            }
                        }
                    }
                }

###############################################################################################
####   Hacky VarScan stuff to get the FET format field included..

                #if FET is defined for file 2, bring it in 
                if($col[8] =~ m/FET/){
                    my @format = split /\:/, $col[8];
                    my $idx=scalar(@format)-1;
                    for my $num (0..(scalar(@format)-1)){
                        if( $format[$idx] =~ m/FET/ ){
                            $idx = $num;
                            last;
                        }
                    }
                    my @values = split /\:/, $col[9];
                    my $fet_value = $values[$idx];
                    @{$varHash{$chr}{$id}}[8] .= ":FET";
                    @{$varHash{$chr}{$id}}[9] .= ":".$fet_value;
                }

##############################################################################################

            } else {

                #add source id
                if(defined($source_ids)){
                    if($col[7] eq "."){
                        $col[7] = "VC=" . $vcfnames[$i];
                    } else {
                        $col[7] = $col[7] . ";VC=" . $vcfnames[$i];
                    }
                }                
                #add to the hash
                @{$varHash{$chr}{$id}} = @col;

                    
            }
            $prev_chr=$chr if $chr;

            
        }
        $self->debug_message("Finished with $prev_chr, printing to file...");
        $self->print_chromosome($prev_chr, $output_file, $varHash{$prev_chr},  \@vcffiles, $passHash{$prev_chr});
    }
    return 1;
}

sub print_chromosome { 
    my $self = shift;
    my $chr = shift;
    my $output_file = shift;
    my $varhash_ref = shift;
    my %varHash = %{$varhash_ref};
   my $vcffiles = shift;
   
   my $passhash_ref = shift; ##only filled for intersect
   my %passHash;
   if($self->require_all_passing && $passhash_ref) {
       %passHash = %{$passhash_ref};
   }   

    #sort by chr, start for clean output
    sub keySort{
        my($x,$y) = @_;
        my @x1 = split(":",$x);
        my @y1 = split(":",$y);
        return($x1[0] <=> $y1[0]);
    }
    my @sortedKeys = sort { keySort($a,$b) } keys %varHash;

    #output
    open(OUTFILE, ">>$output_file") or die "Can't open output file: $!\n";





    #then the body
    foreach my $key (@sortedKeys){

        if($self->require_all_passing){
            #remove lines that aren't passing in all files
            if (@{$varHash{$key}}[6] eq "PASS"){
                if(!(exists($passHash{$key}))){
                    @{$varHash{$key}}[6] = "intersect";
                } else {
                    if ($passHash{$key} < @$vcffiles){
                        @{$varHash{$key}}[6] = "intersect";
                    }
                }
            }
        }        
        print OUTFILE join("\t",@{$varHash{$key}}) . "\n";
    }

    return 1;
}


sub print_header {
    my $self = shift;
    my $output_file = shift;
    my $headers = shift;
    my $newInfo=shift;;
    my $newFilters=shift;
    open(OUTFILE, ">>$output_file") or die "Can't open output file: $!\n";
    
    foreach my $line (@$headers){        
        if ($line =~ /^#CHROM/){
            #dump the info lines just before the column header
            foreach my $line2 (@$newInfo){
                print OUTFILE $line2 . "\n";
            }
            foreach my $line2 (@$newFilters){
                print OUTFILE $line2 . "\n";
            }
        }

        #just output
        print OUTFILE $line . "\n";
    }
}
