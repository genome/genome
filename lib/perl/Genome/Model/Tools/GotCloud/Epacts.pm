package Genome::Model::Tools::GotCloud::Epacts;

use strict;
use warnings;

use Genome;
use Memoize;

class Genome::Model::Tools::GotCloud::Epacts {
    is => ['Genome::Command::Base', 'Genome::Model::Tools::GotCloud'],
    has => [
        vcf => {
            is => 'path',
        },
        ped => {
            is => 'path',
        },
        out => {
            is => 'path',
        },
        type => {
            is => 'text',
            valid_values => [qw(single group make-kin make-group anno)]
        },
        test => {
            is => 'text',
            is_optional => 1,
        },
        pheno => {
            is => 'text',
        },
        cov => {
            is => 'text',
            is_many => 1,
            is_optional => 1,
        },
        min_maf => {
            is => 'integer',
            is_optional => 1,
        },
        max_maf => {
            is => 'integer',
            is_optional => 1,
        },
        min_mac => {
            is => 'integer',
            is_optional => 1,
        },
        min_callrate => {
            is => 'integer',
            is_optional => 1,
        },
        groupf => {
            is => 'text',
            is_optional => 1,
        },
        kin => {
            is => 'text',
            is_optional => 1,
        },
        chr => {
            is => 'text',
            is_many => 1,
            is_optional => 1,
            valid_values => [qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)]
        },
        sepchr => {
            is => 'text',
            is_optional => 1,
        },
        anno => {
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
        field => {
            is => 'text',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    $self->validate;
    my $vcf = $self->vcf;
    my $ped = $self->ped;
    my $out = $self->out;
    my $type = $self->type;
    my $test = $self->test;
    my $pheno = $self->pheno;
    my $epacts_path = $self->epacts_path;
    my $cmd = "$epacts_path $type --vcf $vcf --ped $ped --test $test --pheno $pheno --out $out --run 1";
    if(defined $self->cov){
        $cmd .= join " ", map{" --cov $_ "} $self->cov;
    }
    if(defined $self->min_maf){
        $cmd .=  " --min-maf ".$self->min_maf;
    }
    if(defined $self->max_maf){
        $cmd .= " --max-maf ".$self->max_maf;
    }
    if(defined $self->min_mac){
        $cmd .= " --min-mac ".$self->min_mac;
    }
    if(defined $self->min_callrate){
        $cmd .=" --min-callrate ".$self->min_callrate;
    }
    if(defined $self->chr){
        $cmd .=" --chr ". join " ", $self->chr;
    }
    if(defined $self->sepchr){
        $cmd .=" --sepchr";
    }
    if(defined $self->groupf){
        $cmd .= " --groupf ".$self->groupf;
    }
    if(defined $self->kin){
        $cmd .= " --kinf ".$self->kin;
    }
    if(defined $self->which_skat){
        $cmd .= " --".$self->which_skat;
    }
    if(defined $self->unit){
        $cmd .=" --unit ".$self->unit;
    }
    if(defined $self->anno){
        $cmd .=" --anno";
    }
    if(defined $self->field){
        $cmd .=" --field ".$self->field;
    }
    Genome::Sys->shellcmd(cmd => "$cmd");
    return 1;
}

sub validate {
    my $self = shift;
    Genome::Sys->create_directory($self->out);
    if(defined $self->ped){
        Genome::Sys->validate_file_for_reading($self->ped);
    }
    Genome::Sys->validate_file_for_reading($self->vcf);
    if($self->type eq "group" && !defined $self->groupf){
        die $self->error_message("You must give a group marker file for a group test");
    }
    if($self->type eq "make-kin" && !defined $self->min_maf){
        die $self->error_message("You must give a minimum maf value for make-kin");
    }
    if($self->type eq "make-kin" && !defined $self->min_callrate){
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

