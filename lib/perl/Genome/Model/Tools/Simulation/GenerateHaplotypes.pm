package Genome::Model::Tools::Simulation::GenerateHaplotypes;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
class Genome::Model::Tools::Simulation::GenerateHaplotypes {
    is => 'Command',
    has_optional_input=> [
    output_directory => {
        is =>'Text',
        is_optional=>1,

    },
    reference_panel=>
    {
        is=>'Text',
        is_optional=>0,
        default=>"/gscmnt/sata921/info/medseq/files_for_use_with_simulation_pipeline/hapgen2/CEU.0908.impute.files/CEU.0908.chr22.hap",
    },
    legend_file=>
    { 
        is=>'Text',
        is_optional=>0,
        default=>"/gscmnt/sata921/info/medseq/files_for_use_with_simulation_pipeline/hapgen2/CEU.0908.impute.files/CEU.0908.chr22.legend",
    },
    recombination_map_file=>
    {
        is=>'Text',
        is_optional=>0,
        default=>"/gscmnt/sata921/info/medseq/files_for_use_with_simulation_pipeline/hapgen2/CEU.0908.impute.files/genetic_map_chr22_combined_b36.txt",
    },
    output_file=> {
        is=>'Text',
        is_optional=>0,
        default=>"hapgen2_output",
        is_output=>1,
    },
    number_of_cases => {
        is=>'Number',
        is_optional=>1,
        default=>'1',
    },
    number_of_controls => {
        is=>'Number',
        is_optional=>1,
        default=>'1',
    },
    ],
    has_input => [
    disease_snps=> {
        is=>'Text',
        is_optional=>0,
        doc=>"file of whitespace separated <pos1> <allele1> <rr11> <rr12>  disease allels (pos_in_legend_file, 0 or 1, het risk, hom risk)",
    },
    ],
    has_transient_optional=> [
    _legend_file=> {
        is_output=>1,
        is=>'Text',
    },
    _cases_file=> {
        is_output=>1,
        is=>'Text',
    },
    _controls_file=> {
        is_output=>1,
        is=>'Text',
    },
    ],
};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}

sub execute {
    my $self = shift;
    my $disease_allele_string = $self->process_disease_snps_file($self->disease_snps);
    unless($disease_allele_string) {
        $self->error_message("disease snps file could not be properly parsed!");
        return 0;
    }
    my( $ref_panel, $legend_file, $recomb_file, $output_file, $cases, $controls) = ($self->reference_panel, $self->legend_file, $self->recombination_map_file, $self->output_file, $self->number_of_cases, $self->number_of_controls);
    if($self->output_directory) {
        $output_file = $self->output_directory . "/$output_file";
    }
    my $cmd = "/gscmnt/sata921/info/medseq/files_for_use_with_simulation_pipeline/hapgen2/hapgen2 ";
    $cmd .= " -h $ref_panel ";
    $cmd .= " -l $legend_file ";
    $cmd .= " -m $recomb_file ";
    $cmd .= " -o $output_file "; 
    $cmd .= "-n $controls $cases ";
    $cmd .= " $disease_allele_string ";
    print "trying to execute $cmd \n";
    my $rv = Genome::Sys->shellcmd(cmd=>"$cmd");
    $self->_legend_file($output_file .".legend");
    $self->_cases_file($output_file . ".cases.haps");
    $self->_controls_file($output_file . ".controls.haps");
    return $rv;
}



sub process_disease_snps_file {
    my $self = shift;
    my $disease_snps_file = shift;
    my $fh = Genome::Sys->open_file_for_reading($disease_snps_file);
    my $string = "-dl ";
    while(my $line = $fh->getline) {
        chomp($line);
        $string .= " $line ";
    }
    return $string;
}
1;
