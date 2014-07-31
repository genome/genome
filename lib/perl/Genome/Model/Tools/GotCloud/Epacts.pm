package Genome::Model::Tools::GotCloud::Epacts;

use strict;
use warnings;

use Genome;
use Memoize;

class Genome::Model::Tools::GotCloud::Umake {
    is => ['Genome::Command::Base', 'Genome::Model::Tools::GotCloud'],
    has => [
        vcf_file => {
            is => 'path',
        },
        pedigree_file => {
            is => 'path',
        },
        output_directory => {
            is => 'path',
        },
        type => {
            is => 'text',
            valid_values => [qw (single group make-kin make-group anno)]
        },
        test => {
            is => 'text',
            is_optional => 1,
        },
        phenotype => {
            is => 'text',
        },
        covariates => {
            is => 'text',
            is_many => 1,
        },
        minimum_maf => {
            is => 'integer',
            is_optional => 1,
        },
        maximum_maf => {
            is => 'integer',
            is_optional => 1,
        },
        minimum_mac => {
            is => 'integer',
            is_optional => 1,
        },
        minimum_callrate => {
            is => 'integer',
            is_optional => 1,
        },
        marker_group_file => {
            is => 'text',
            is_optional => 1,
        },
        kinship_matrix => {
            is => 'text',
            is_optional => 1,
        },
        chromosomes => {
            is => 'text',
            is_many => 1,
            is_optional => 1,
            valid_values => [qw (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)]
        },
        separate_by_chr => {
            is => 'text',
            is_optional => 1,
        },
        annotate => {
            is => 'text',
            is_optional => 1,
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            is_optional => 1,
        },
        which_skat => {
            is => 'text',
            is_optional => 1,
        },
        unit => {
            is => 'integer',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    $self->validate;
    my $vcf = $self->vcf_file;
    my $ped = $self->pedigree_file;
    my $out = $self->output_directory;
    my $type = $self->type;
    my $test = $self->test;
    my $pheno = $self->pheno;
    my $cov = $self->covariates;
    my $min_maf = $self->minimum_maf;
    #my $max_maf = $self->maximum_maf;
    #my $min_mac = $self->minimum_mac;
    my $min_call = $self->minimum_callrate;
    #my $grp = $self->marker_group_file;
    #my $kin = $self->kinship_matrix;
    my $chr = $self->chromosomes;
    my $sep_by_chr = $self->separate_by_chr;
    my $anno = $self->annotate;
    my $ref = $self->reference_build;
    #my $skat_test = $self->which_skat;
    my $epacts_path = $self->epacts_path;
    my $unit = $self->unit;
    #Genome::Sys->shellcmd(cmd =>"$epacts_path");
    if($type eq "make-kin"){
        my $cmd;
        if(defined $self->chr){
            $cmd = "$epacts_path $type --vcf $vcf --ped $ped --min-maf $min_maf --min-callrate $min_call --out $out --chr $chr";
        }
        if(defined $self->separate_by_chr){
            $cmd = "$epacts_path $type --vcf $vcf --ped $ped --min-maf $min_maf --min-callrate $min_call --out $out --sepchr";
        }
        else{
            $cmd = "$epacts_path $type --vcf $vcf --ped $ped --min-maf $min_maf --min-callrate $min_call --out $out";
        }
        Genome::Sys->shellcmd(cmd => "$cmd");
    }
    my $cmd = "$epacts_path $type --vcf $vcf --ped $ped --test $test --pheno $pheno --out $out";
    if(defined $self->covariates){
        $cmd .= join " ", map{"--cov $_"} $self->covariates;
    }
    if(defined $self->minimum_maf){
        $cmd .=  "--min-maf".$self->minimum_maf;
    }
    if(defined $self->maximum_maf){
        $cmd .= "--max-maf".$self->maximum_maf;
    }
    if(defined $self->minimum_callrate){
        $cmd .="--min-callrate".$self->minimum_callrate;
    }
    if(defined $self->marker_group_file){
        $cmd .= "--groupf".$self->marker_group_file;
    }
    if(defined $self->kinship_matrix){
        $cmd .= "--kinf".$self->kinship_matrix;
    }
    if(defined $self->skat_test){
        $cmd .= "--".$self->skat_test;
    }
    if(defined $self->unit){
        $cmd .="--".$self->unit;
    }
    Genome::Sys->shellcmd(cmd => "$cmd");
    return 1;
}

sub validate {
    my $self = shift;
    Genome::Sys->create_directory($self->output_directory);
    if(defined $self->pedigree_file){
        Genome::Sys->validate_file_for_reading($self->pedigree_file);
    }
    Genome::Sys->validate_file_for_reading($self->vcf_file);
    if($self->type eq "group" && !defined $self->marker_group_file){
        die $self->error_message("You must give a group marker file for a group test");
    }
    if($self->type eq "make-kin" && !defined $self->min_maf){
        die $self->error_message("You must give a minimum maf value for make-kin");
    }
    if($self->type eq "make-kin" && !defined $self->minimum_callrate){
        die $self->error_message("You must give a minimum call rate value for make-kin");
    }
    if(defined $self->anno && !defined $self->reference_build){
        die $self->error_message("You must give the reference build to annotate");
    }
    if($self->test eq "skat" && !defined $self->which_skat){
        die $self->error_message("You must specify WHICH skat test you want to run if you choose skat");
    }
    return 1;
}

