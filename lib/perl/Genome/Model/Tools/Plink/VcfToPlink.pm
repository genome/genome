package Genome::Model::Tools::Plink::VcfToPlink;     # rename this when you give the module file a different name <--


use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use IO::File;
use Genome::Sys;

class Genome::Model::Tools::Plink::VcfToPlink {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
    vcf_file	=> { is => 'Text', doc => "Input VCF file to be converted" , is_optional => 0, is_input => '1'},
    pedigree_file	=> { is => 'Text', doc => "A standard pedigree file with phenotype, gender, & relationships" , is_optional => 0, is_input => '1'},
    phenotype_file	=> { is => 'Text', doc => "Name for output phenotype file for PLINK" , is_optional => 0, is_input => '1'},
    plink_file	=> { is => 'Text', doc => "Name for output PLINK file" , is_optional => 0, is_input => '1'},
    output_directory	=> { is => 'Text', doc => "Name of output directory", is_optional => 0, is_input => '1'},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Converts a VCF to PLINK format"                 
}

sub help_synopsis {
    return <<EOS
This command converts a VCF to PLINK format
EXAMPLE:	gmt plink vcf-to-plink --vcf-file [file] --pedigree-file [file] --phenotype-file [name] --plink-file [name] --output-directory [name]...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
The steps to convert a VCF to PLINK are as follows:
1.) Run VCFtools to convert to PLINK BED format
2.) Update individual & family IDs
3.) Update genders
4.) Update parent relationships
5.) Create a correct phenotype file to be used with --pheno command in PLINK
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my ($self) = @_;

    my $input_vcf = $self->vcf_file;
    my $input_ped = $self->pedigree_file;
    my $output_pheno = $self->phenotype_file;
    my $output_plink = $self->plink_file;
    my $output_dir = $self->output_directory;
    my $convert_cmd = "/gscmnt/ams1161/info/model_data/kmeltzst/Software/vcftools_0.1.9/bin/vcftools --gzvcf $input_vcf --plink --out $output_dir/tmp";
    Genome::Sys->shellcmd(cmd=>$convert_cmd);
    
    my $out_family_file = IO::File->new("$output_dir/family_file.txt", ">");
    my $out_sex_file = IO::File->new("$output_dir/sex_file.txt", ">");
    my $out_parents_file = IO::File->new("$output_dir/parents_file.txt", ">");
    my $out_pheno_file = IO::File->new($output_dir."/".$output_pheno, ">");
    my $in_file = IO::File->new($input_ped);
    while(my $line = $in_file->getline){
        chomp($line);
        my ($family_id,$ind_id,$father_id,$mother_id,$sex,$pheno)=(split/\t/,$line);
        $out_family_file->print("$ind_id\t$ind_id\t$family_id\t$ind_id\n");
        $out_sex_file->print("$family_id\t$ind_id\t$sex\n");
        $out_parents_file->print("$family_id\t$ind_id\t$father_id\t$mother_id\n");
        $out_pheno_file->print("$family_id\t$ind_id\t$pheno\n");
    }

    my $update_id_cmd = Genome::Model::Tools::Plink->path_to_binary() . " --file $output_dir/tmp --update-ids $output_dir/family_file.txt --make-bed --out $output_dir/tmp2";
    Genome::Sys->shellcmd(cmd=>$update_id_cmd);

    my $update_parents_cmd = Genome::Model::Tools::Plink->path_to_binary() . " --bfile $output_dir/tmp2 --update-parents $output_dir/parents_file.txt --make-bed --out $output_dir/tmp3";
    Genome::Sys->shellcmd(cmd=>$update_parents_cmd);

    my $update_sex_cmd = Genome::Model::Tools::Plink->path_to_binary() . "--bfile $output_dir/tmp3 --update-sex $output_dir/sex_file.txt --make-bed --out $output_dir/$output_plink";
    Genome::Sys->shellcmd(cmd=>$update_sex_cmd);

    my $cleanup_cmd = "rm -f tmp*";
    #Genome::Sys->shellcmd(cmd=>$cleanup_cmd);

    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;
