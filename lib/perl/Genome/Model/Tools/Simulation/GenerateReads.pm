package Genome::Model::Tools::Simulation::GenerateReads;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
class Genome::Model::Tools::Simulation::GenerateReads {
    is => 'Command',
    has_optional_input => [
    diploid_fasta => {
        type => 'String',
        is_optional => 0,
        is_input=>1,
        doc=>"you can get this from alter-reference-sequences if you don't know wtf you're doing and you're in this dir",
    },
    mason_parameters => {
        type => 'String',
        is_optional => 0,
        doc => 'string of parameters to pass to mason',
        default=> " -i -sq -n 100 -mp -ll 500 ",
    },
    output_bam_name => { 
        type => 'String',
        doc => 'the name of the output bam that will contain all your simulated reads',
        is_optional=>1,
    },
    rename_reads=> {
        is_optional=>1,
        default=>1,
        doc=>'rename the reads to be more useful for downstream mapping/other analysis',
    },
    coverage_level=> {
        is_optional=>1,
        default=>100,
        doc=>'take the supplied fasta length and (assuming 100bp reads) generate the required number of reads to get to this haploid coverage number',
    },
    ],
};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}

sub execute {
    my $self=shift;
    #FIXME do this as a calc maybe
    my $output_bam = $self->diploid_fasta;
    $output_bam =~s/\.fasta//g;
    $self->output_bam_name($output_bam);
    ####

    $DB::single=1;
    my $mason_cmd="/gscuser/charris/bin/mason illumina " . $self->mason_parameters;  #FIXME deploy, no hardcode
    my $current_dir = getcwd;
    
    my $temp_dir = Genome::Sys->base_temp_directory();
    chdir($temp_dir);
    $mason_cmd.= " -o mason_output ";
    $mason_cmd.=$self->diploid_fasta;
    my $num_reads_to_generate = $self->estimate_reads_for_coverage_target($self->coverage_level, $self->diploid_fasta, 100);
    $mason_cmd.=" -N $num_reads_to_generate ";

    #FIXME: add in the coverage parameter to the mason command after calculating it
    #FIXME: another fixme just to draw attention to this glaring omission
    unless(Genome::Sys->shellcmd(cmd=>$mason_cmd)) {
        $self->error_message("Error running mason!");
        return 0;
    }
    my $fixed_sam_file = $self->imbue_read_names_with_more_meaning("mason_output.sam");
    my $final_bam = $self->output_bam_name;
    $final_bam =~ s/.bam//; #samtools adds this like a mofo so take it out to make the user think he is doing something rather than try to epxlain how dumb he is to himself
    chdir($current_dir);
    my $make_my_bam = "samtools view $temp_dir/$fixed_sam_file -u -S | samtools sort -n - $final_bam";  #FIXME No samtools hardcode;
    unless(Genome::Sys->shellcmd(cmd=>$make_my_bam)) {
        $self->error_message("Error turning improved sam into bam, get heng li on the batphone.");
        return 0;
    }

    

}
    1;

sub imbue_read_names_with_more_meaning {
    my $self=shift;
    my $sam_file = shift;
    return $sam_file; #FIXME: do something more than nothing in this method
}

sub estimate_reads_for_coverage_target {
    my $self = shift;
    my $desired_coverage = shift;
    my $diploid_fasta = shift;
    my $read_length = shift;
    my @fasta_headers = `grep "^>" $diploid_fasta`;  
    my $total_length;
    for my $header (@fasta_headers) {
        my($fast_name, $fasta_size, $desc) = split "\t", $header;
        $total_length+=$fasta_size;
    }
    $total_length/=2; #we made two haplotypes of every region so divide our summed length by 2 to correct
    my $num_reads_to_generate = ceil(($total_length*$desired_coverage)/($read_length*2));
    return $num_reads_to_generate;
}




