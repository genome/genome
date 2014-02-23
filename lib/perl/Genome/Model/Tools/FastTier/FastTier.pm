package Genome::Model::Tools::FastTier::FastTier;

use strict;
use warnings;

use Genome;     
use File::Basename;
use File::Copy;

class Genome::Model::Tools::FastTier::FastTier {
    is => 'Command',
    has => [
        variant_bed_file => {
            type => 'Text',
            is_input => 1,
        },
        tier_file_location => {
            type => 'Text',
            is_input => 1,
            doc => 'Use this to point to a directory containing tier1.bed - tier4.bed in order to use different bed files for tiering',
        },
        skip_line_count => {
            type => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'Set this if the input has duplicates',
        },
    ],
    has_optional => [
        tiering_version => {
            is => 'Text',
            doc => "number specifying the version of bed files to use",
            is_input => 1,
            default => 3,
        },
        tier1_output => {
            is => 'Text',
            is_optional => 1,
            is_input => 1, 
            is_output => 1,
        },
        tier2_output => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        tier3_output => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        tier4_output => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        intersect_bed_bin_location => {
            type => 'Text',
            is_input => 1,
            default => '/gsc/bin/intersectBed-pipes',  #$ENV{GENOME_SW} . '/bedtools/installed-64/intersectBed',
            doc => 'The path and filename of intersectBed',
        },
        output_directory => {
            is => 'Text',
            doc => 'Set this to drop output files somewhere other than the location of the input',
        },
        _tier1_bed => {
            type => 'Text',
        },
        _tier2_bed => {
            type => 'Text',
        },
        _tier3_bed => {
            type => 'Text',
        },
        _tier4_bed => {
            type => 'Text',
        }, 
    ]

};

sub sub_command_sort_position { 15 }

sub help_brief {
    "This tool uses the GC customized version of intersectBed in order to stream a bedfile end to end, dropping the appropriate items in the appropriate tier files."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt fast-tier fast-tier
EOS
}

sub execute {
    my $self=shift;
    unless(-s $self->variant_bed_file) {
        $self->error_message("The variant file you supplied: " . $self->variant_bed_file . " appears to be 0 size. You need computer more better.");
        return;
    }
    if($self->output_directory){
        my $output_directory = $self->output_directory;
        unless(-d $output_directory){
            die $self->error_message("Could not locate output directory!");
        }
        my $base = basename($self->variant_bed_file);
        for my $tier(1..4){
            my $output_accessor = "tier".$tier."_output";
            unless ($self->$output_accessor){
                $self->$output_accessor($output_directory."/".$base.".tier".$tier);
            }
        }
    } else {
        for my $tier(1..4){
            my $output_accessor = "tier".$tier."_output";
            unless ($self->$output_accessor){
                $self->$output_accessor($self->variant_bed_file.".tier".$tier);
            }
        }
    }

    #if the user specified an alternate location for tier bed files, check and load them
    if(defined($self->tier_file_location)){
        unless(-d $self->tier_file_location){
            $self->error_message("You must specify a directory containing the 4 tier files.");
            die $self->error_message;
        }
        my @tiers = map { $self->tier_file_location."/tier".$_.".bed";} (1,2,3,4);
        for my $t (@tiers){
            # This previously checked -s, but some things (like mouse) do not have every tier available
            unless(-e $t){
                $self->error_message("Could not locate a bed file at ".$t."\n");
                die $self->error_message;
            }
        }
        $self->_tier1_bed($tiers[0]);
        $self->_tier2_bed($tiers[1]);
        $self->_tier3_bed($tiers[2]);
        $self->_tier4_bed($tiers[3]);
    }

    $self->debug_message("Using tier 1 bed file  at ".$self->_tier1_bed ."\n");
    $self->debug_message("Using tier 2 bed file  at ".$self->_tier2_bed ."\n");
    $self->debug_message("Using tier 3 bed file  at ".$self->_tier3_bed ."\n");
    $self->debug_message("Using tier 4 bed file  at ".$self->_tier4_bed ."\n");

    if($self->tiering_version >= 3){

        my $sorted_input = Genome::Sys->create_temp_file_path;
        my $sort_incoming = Genome::Model::Tools::Joinx::Sort->create(input_files => [$self->variant_bed_file],output_file => $sorted_input);
        unless($sort_incoming->execute){
            die $self->error_message("Could not call sort on input: ".$self->variant_bed_file);
        }
        my $insertions = Genome::Sys->create_temp_file_path;

        my @tier_temp_files;
        map{ push @tier_temp_files, Genome::Sys->create_temp_file_path;}(1..4);

        my $tier1_cmd = "joinx intersect -f -o ".$tier_temp_files[0]." -a ".$sorted_input." -b ".$self->_tier1_bed." --miss-a - ";
        my $tier2_cmd = "joinx intersect -f -o ".$tier_temp_files[1]." -a - -b ".$self->_tier2_bed." --miss-a - ";
        my $tier3_cmd = "joinx intersect -f -o ".$tier_temp_files[2]." -a - -b ".$self->_tier3_bed." --miss-a - ";
        my $tier4_cmd = "joinx intersect -f -o ".$tier_temp_files[3]." -a - -b ".$self->_tier4_bed." --miss-a - "; 

        my $insertions_cmd = $insertions;

        my $cmd = $tier1_cmd . " | " . $tier2_cmd . " | " . $tier3_cmd . " | " . $tier4_cmd . " > " . $insertions_cmd;

        unless( Genome::Sys->shellcmd( cmd => $cmd ) ){
            die $self->error_message("Fast tier command did not complete!");
        }

        my @ins_temp_files;
        map{ push @ins_temp_files, Genome::Sys->create_temp_file_path;}(1..4);

        if(-s $insertions){

            my $tier4_ins_cmd = "joinx intersect --adjacent-insertions -f -o ".$ins_temp_files[3]." -a ".$insertions." -b ".$self->_tier4_bed." --miss-a -";
            my $tier3_ins_cmd = "joinx intersect --adjacent-insertions -f -o ".$ins_temp_files[2]." -a - -b ".$self->_tier3_bed." --miss-a - ";
            my $tier2_ins_cmd = "joinx intersect --adjacent-insertions -f -o ".$ins_temp_files[1]." -a - -b ".$self->_tier2_bed." --miss-a - ";
            my $tier1_ins_cmd = $ins_temp_files[0];

            my $ins_cmd = $tier4_ins_cmd . " | " . $tier3_ins_cmd . " | " . $tier2_ins_cmd . " > " . $tier1_ins_cmd;

            unless( Genome::Sys->shellcmd( cmd => $ins_cmd ) ){
                die $self->error_message("Fast tier command for insertions did not complete!");
            }


            for my $tier_num (1..4){
                my $output = "tier".$tier_num."_output";
                $tier_num--;
                my $sort = Genome::Model::Tools::Joinx::Sort->create(input_files => [$tier_temp_files[$tier_num],$ins_temp_files[$tier_num]],output_file => $self->$output, merge_only => 1);
                unless($sort->execute){
                    die $self->error_message("Could not merge insertions with the rest of the results.");
                }
            }
        } 
        else {

            for my $tier_num(1..4){
                my $output = "tier".$tier_num."_output";
                $tier_num--;
                Genome::Sys->copy_file($tier_temp_files[$tier_num], $self->$output);
            }
        }
    } 
    else {
        my $tier1_cmd = $self->intersect_bed_bin_location." -wa -vf stdout -of ".$self->tier1_output." -u -a ".$self->variant_bed_file." -b ".$self->_tier1_bed;
        my $tier2_cmd = $self->intersect_bed_bin_location." -wa -vf stdout -of ".$self->tier2_output." -u -a stdin -b ".$self->_tier2_bed;
        my $tier3_cmd = $self->intersect_bed_bin_location." -wa -vf stdout -of ".$self->tier3_output." -u -a stdin -b ".$self->_tier3_bed;
        my $tier4_cmd = $self->tier4_output;

        my $cmd = $tier1_cmd . " | " . $tier2_cmd . " | " . $tier3_cmd . " > " . $tier4_cmd;

        unless( Genome::Sys->shellcmd( cmd => $cmd ) ){
            die $self->error_message("Fast tier command did not complete!");
        }
    }
    unless($self->skip_line_count){
        my $input = $self->line_count($self->variant_bed_file);
        my $output=0;
        for my $tier (1..4){
            my $out = "tier".$tier."_output";
            $output += $self->line_count($self->$out);
        }
        unless(($input - $output)==0){
            die $self->error_message("Lines of input: ".$input."\t did not match lines of output: ".$output."\n");
        }
    }
    return 1;
}


sub line_count {
    my $self = shift;
    my $input = shift;
    unless( -e $input ) {
        die $self->error_message("Could not locate file for line count: $input");
    }
    my $result = `joinx sort -u $input | wc -l`; 
    my ($answer)  = split /\s/,$result;
    return $answer
}
