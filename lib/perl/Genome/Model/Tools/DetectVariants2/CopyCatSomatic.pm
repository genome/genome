package Genome::Model::Tools::DetectVariants2::CopyCatSomatic;

use strict;
use warnings;

use Cwd;
use Genome;
use Workflow::Simple;

class Genome::Model::Tools::DetectVariants2::CopyCatSomatic{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
        tumor_window_dir => {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'dir containing tumor window file',
        },
        normal_window_dir => {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'dir containing normal window file',
        },
        per_library => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            default => 1,
            doc => 'do per-library correction',
        },
        per_readlength => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            default => 1,
            doc => 'do per-readlength correction',
        },
        tumor_samtools_file => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            default => "NULL",
            doc => 'path to tumor samtools file',
        },
        reference_build_id => {
            type => 'String',
            is_optional => 0,
            doc => 'reference build id',
        },
        annotation_directory => {
            type => 'Path',
            is_optional => 0,
            is_input => 1,
            doc => 'annotation data',
        },
        lsf_resource => {
            default_value => "-R 'rusage[mem=4000] select[type==LINUX64 && maxtmp>10000] span[hosts=1]' -M 4000000 -n 4",
        },
     ],
};


sub _detect_variants {
    my $self = shift;
    
    my $cmd = Genome::Model::Tools::CopyCat::Somatic->create( 
        tumor_window_file => $self->tumor_window_dir . "/readcounts.wind",
        normal_window_file => $self->normal_window_dir . "/readcounts.wind",
        output_directory => $self->_temp_staging_directory,
        per_library => $self->per_library,
        per_read_length => $self->per_readlength,
        processors => 4,
        genome_build => $self->reference_name,
        tumor_samtools_file => $self->tumor_samtools_file,
        annotation_directory => $self->annotation_directory,
        );

    unless($cmd->execute){
        $self->error_message("Failed to run CopyCat command.");
        die $self->error_message;
    }
    
    my $cnvs = $self->_temp_staging_directory."/cnvs.hq";    
    my $outfile = $self->_temp_staging_directory."/segs.paired.dat";    
    system("ln -s $outfile $cnvs");

    return 1;
}



sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

sub _sort_detector_output {
    return 1;
}

sub reference_name {
    my $self = shift;

    #get reference genome name
    my $refbuild_id = $self->reference_build_id;
    unless($refbuild_id){
        die $self->error_message("Received no reference build id.");
    }
    print "refbuild_id = ".$refbuild_id."\n";
    my $ref_seq_build = Genome::Model::Build->get($refbuild_id);
    unless($ref_seq_build){
        die $self->error_message("Not a valid reference build id");
    }
    my $ref_seq_name = $ref_seq_build->name;
    
    if ($ref_seq_name eq "NCBI-human-build36"){
        return("36");
    } elsif ($ref_seq_name eq "GRCh37-lite-build37"){
        return("37");
    } elsif ($ref_seq_name eq "UCSC-mouse-buildmm9"){
        return("mm9");
    } else {
        die $self->error_message("this tool only supports NCBI-human-build36, GRCh37-lite-build37, and UCSC-mouse-buildmm9");
    }
}

1;
