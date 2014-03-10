package Genome::Model::Tools::Snp::SomaticCaptureDbsnp;

use strict;
use warnings;

use Genome;
use Command;

class Genome::Model::Tools::Snp::SomaticCaptureDbsnp{
    is => 'Command',
    has => [
        dbsnp_snvs_file => {
            type => 'String',
            doc => 'This is the file to contain snvs that are in dbSNP',
        },
        novel_snvs_file => {
            type => 'String',
            doc => 'This is the file to contain snvs that are NOT in dbSNP',
        },
        snv_file => { 
            type => 'String',
            doc => 'snv file',
        },
        dbsnp_build_id => {
            type => 'String',
            is_optional => 1,
            doc => "Build id of the dbSNP build to intersect with",
        },
        dbsnp_bed_file => {
            type => 'String',
            is_optional => 1,
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
    unless(-e $self->snv_file) {
        $self->error_message("Snps file is not a file: " . $self->snv_file);
        return;
    }
    my $snv_bed_path = $self->dbsnp_bed_file;
    my $snv_input_path = $self->convert_to_bed($self->snv_file);
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

    print "snv_input_path = ".$snv_input_path."\n";
    print "snv_bed_path = ".$snv_bed_path."\n";
    print "novel_snvs_file = ".$novel_snvs_file."\n";
    print "dbsnp_snvs_file = ".$dbsnp_snvs_file."\n";


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

# return a path to a temp file containing a bed version of the snv_file
sub convert_to_bed {
    my $self = shift;
    my $snv_file = shift;
    $self->debug_message("Converting input to dbsnp");
    my $temp_bed_file = Genome::Sys->create_temp_file_path;
    my $temp_sorted_file = Genome::Sys->create_temp_file_path;
    my $temp_headless_file = Genome::Sys->create_temp_file_path;
    my $snv_fh = Genome::Sys->open_file_for_reading($snv_file);
    my $snv_headless_fh = Genome::Sys->open_file_for_writing($temp_headless_file);
    $snv_fh->getline;
    $snv_fh->getline;
    $snv_fh->getline;
    while(my $line = $snv_fh->getline){
        print $snv_headless_fh $line;
        print $line;
    }
    $snv_fh->close;
    $snv_headless_fh->close;
    my @snv_files = ($temp_headless_file);
    my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create( input_files => \@snv_files, output_file => $temp_sorted_file);

    unless($sort_cmd->execute){
        die $self->error_message("Could not execute joinx sort command.");
    }


    my $snv_anno_fh = Genome::Sys->open_file_for_reading($temp_sorted_file);
    my $temp_fh = Genome::Sys->open_file_for_writing($temp_bed_file);
    while( my $line = $snv_anno_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop,$ref,$var) = split /\s+/, $line;
        print $temp_fh join("\t",($chr,$start-1,$stop,$ref,$var))."\n";
    }
    $temp_fh->close;
    $snv_anno_fh->close;
    return $temp_bed_file;
}

1;
