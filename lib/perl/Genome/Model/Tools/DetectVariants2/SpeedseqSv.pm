#SpeedseqSv.pm


package Genome::Model::Tools::DetectVariants2::SpeedseqSv;
use warnings;
use strict;

use Genome;
use File::Basename;
use IPC::System::Simple;
use Genome::File::Tsv;

class Genome::Model::Tools::DetectVariants2::SpeedseqSv {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has_optional => [
        _legend_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ File::Spec->join($_temp_staging_directory, 'legend.tsv'); },
        },
    ]
};

	my $temp_directory = Genome::Sys->create_temp_directory();
	my $pkg = 'Genome::Model::Tools::Speedseq::Sv';

# For the Parameters what I need to do is make two hashes.
# The first one will have the name of the command in the SV class and the letter calling it.  
# The second Hash will have the class string name and the value wil be either a boolean or a file location.
	my %ref_library = (
		B => 'full_bam_file',
		S => 'split_read_bam_file',
		D => 'discordant_read_bam_file',
		R => 'reference_fasta',
		g => 'genotype_svtyper',
		d => 'CNVnator_read_depth',
		A => 'annotate_vcf',
		k => 'keep_temp_files',
		#End of Boolean files#
		t => 'threads',
		x => 'exclude_bed_file',
		m => 'min_sample_weight',
		r => 'trim_threshold',
		T => 'temp_directory',	
		K => 'config_file',
		v => 'verbose',
	);


sub _detect_variants {
	my $self = shift;
	
	my @fullBam = ($self->aligned_reads_input,$self->control_aligned_reads_input);

	my $aligned_bams = join(',',@fullBam);

	my %final_cmd = (
   		output_prefix => $self->_sv_staging_output,
   		full_bam_file => $aligned_bams,
		$self->find_split(@fullBam),
		$self->find_discord(@fullBam),
	);	

	my %list_params = $self->split_params_to_letter();
 	
	while (my ($key, $value) = each(%ref_library)) {$final_cmd{$value} = $list_params{$key} if exists $list_params{$key};}

	my $set = $pkg->create(%final_cmd);
	$set->execute();
};

sub find_split{
        my $self = shift;
        my @bam_dir = @_;
	my @final = ();
        
	foreach(@bam_dir){
		my $split;
        	my @dir_split = split('/',$_);
        	my $editor = pop(@dir_split);
        	$editor =~ s/aligned/splitters/g;
        	push (@dir_split, $editor);
        	$split = join('/', @dir_split);
        	push (@final, $split);
	}
	my $combined_splits = join (',',@final);
        return (
		split_read_bam_file => $combined_splits,
	);
};

sub find_discord{
        my $self = shift;
	my @bam_dir = @_;
	my @final = ();
	
	foreach (@bam_dir){
		my $discord;		
		my @dir_split = split('/',$_);
		my $editor = pop(@dir_split);
		$editor =~ s/aligned/discordants/g;
		push (@dir_split, $editor);
		$discord = join('/', @dir_split);
		push (@final, $discord);
	}
	my $combined_splits = join (',',@final);
        return (
                discordant_read_bam_file => $combined_splits,
        );

};

sub split_params_to_letter {
        my $self = shift;
        my $parms = $self->params;
        my %params_hash = ();
	my @params = split(',',$parms);
	foreach (@params){
		if ($_ =~ /:/){
			my $num = substr($_,1,1);
			my $value = substr($_,3);
	                $params_hash{$num} = $value; 
		}
                else {
			my $num = substr($_,1,1);
	                $params_hash{$num} = 'true'; 
                }
	}
	return %params_hash;

};


