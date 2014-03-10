package Genome::Model::Tools::Snp::DbsnpIntersection;

use strict;
use warnings;

use Genome;
use Command;

class Genome::Model::Tools::Snp::DbsnpIntersection{
    is => 'Command',
    has => [
        dbsnp_snvs_file => {
            type => 'String',
            is_output => 1,
            doc => 'This is the file to contain snvs that are in dbSNP',
        },
        novel_snvs_file => {
            type => 'String',
            is_output => 1,
            doc => 'This is the file to contain snvs that are NOT in dbSNP',
        },
        snv_bed_file => { 
            type => 'String',
            is_input => 1,
            doc => 'snv file',
        },
        dbsnp_build_id => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => "Build id of the dbSNP build to intersect with",
        },
        dbsnp_bed_file => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => "Specify either a dbsnp bed file OR a dbsnp_build_id, NOT both.",
        },
    ]
};

sub help_detail {
    "This script performs a comparison of a maq cns2snp output file with a dbSNP file. The comparisons are made by position only as a dbSNP file does not include allele information at this time."
}

sub execute {
    my $self=shift;
    unless(defined($self->dbsnp_bed_file) xor defined($self->dbsnp_build_id)){
        die $self->error_message("You must define either a dbsnp bed file or a dbsnp build id, but not both.");
    }
    #Check on the file names
    unless(-e $self->snv_bed_file) {
        $self->error_message("Snps file is not a file: " . $self->snv_bed_file);
        return;
    }
    my $snv_bed_path = $self->dbsnp_bed_file;
    my $snv_input_path = $self->sort_bed($self->snv_bed_file);
    if(defined($self->dbsnp_build_id)){
        my $dbsnp_build_id = $self->dbsnp_build_id;
        my $dbsnp_build = Genome::Model::Build->get($dbsnp_build_id);
        $snv_bed_path = $dbsnp_build->snvs_bed;
    }
    unless(defined($snv_bed_path)){
        die $self->error_message("Could not locate a path to dbsnp bed file.");
    }
    my $novel_snvs_file = $self->novel_snvs_file;
    my $dbsnp_snvs_file = $self->dbsnp_snvs_file;

    my $snv_compare = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $snv_input_path,
        input_file_b => $snv_bed_path,
        miss_a_file => $novel_snvs_file,
        output_file => $dbsnp_snvs_file,
        exact_pos => 1,
        iub_match => 1,
    );

    unless ($snv_compare->execute){
        die $self->error_message("Couldn't create snv comparison tool!");
    }


    return 1;
}

# return a path to a temp file containing a sorted snv_file
sub sort_bed {
    my $self = shift;
    my $snv_file = shift;
    $self->debug_message("Sorting input snv file");
    my $temp_sorted_file = Genome::Sys->create_temp_file_path;
    my @snv_files = ($snv_file);
    my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create( input_files => \@snv_files, output_file => $temp_sorted_file);

    unless($sort_cmd->execute){
        die $self->error_message("Could not execute joinx sort command.");
    }

    return $temp_sorted_file;
}

1;
