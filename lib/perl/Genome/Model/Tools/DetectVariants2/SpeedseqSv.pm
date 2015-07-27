#SpeedseqSv.pm


package Genome::Model::Tools::DetectVariants2::SpeedseqSv;
use warnings;
use strict;

use Genome;
use File::Basename qw(basename dirname);
use File::Spec;
use IPC::System::Simple;
use Genome::File::Tsv;

class Genome::Model::Tools::DetectVariants2::SpeedseqSv {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
};

	my $temp_directory = Genome::Sys->create_temp_directory();
	my $pkg = 'Genome::Model::Tools::Speedseq::Sv';

# For the Parameters what I need to do is make two hashes.
# The first one will have the name of the command in the SV class and the letter calling it.  
# The second Hash will have the class string name and the value wil be either a boolean or a file location.

sub _detect_variants {
	my $self = shift;
	
	my @fullBam = ($self->aligned_reads_input,$self->control_aligned_reads_input);

	my $aligned_bams = join(',',@fullBam);

	my %library = $self->get_parameter_hash();

	while (my ($key, $value) = each(%library)){
		print "$key => $value\n";
	}
	my %final_cmd = ();

	my %list_params = $self->split_params_to_letter();
 	
	while (my ($key, $value) = each(%library)) {$final_cmd{$value} = $list_params{$key} if exists $list_params{$key};}

	%final_cmd = (
		%final_cmd,
   		output_prefix => $self->_sv_staging_output,
   		full_bam_file => $aligned_bams,
		temp_directory => Genome::Sys->create_temp_directory(),
		$self->find_file("splitters","split_read_bam_file",@fullBam),
		$self->find_file("discordants","discordant_read_bam_file",@fullBam),
	);	

	my $set = $pkg->create(%final_cmd);
	$set->execute();
};

sub find_file{
        my $self = shift;
	my $disc = shift;
	my $value = shift;
	my @bam_dir = @_;
	my @final = ();
	
	foreach (@bam_dir){
		my $editor = basename($_);
		my $dir = dirname($_);
		$editor =~ s/\.bam/\.$disc\.bam/g;
		my $discord = "$dir/$editor";
		if (!-s $discord) {die $self->error_message("File couldn't be found: $discord  Bam Files Must be aligned with Speedseq.")};
		push (@final, $discord);
	}
	my $combined_splits = join (',',@final);
        return (
                $value => $combined_splits,
        );
};

sub split_params_to_letter {
        my $self = shift;
        my $parms = $self->params;
        my %params_hash = ();
	my @params = split(',',$parms);
	foreach (@params){
		if ($_ =~ /:/){
			my @options = split(':',$_);
			my $value = pop(@options);
			my $num = shift(@options);
	                $num =~ s/-//gi;
			$params_hash{$num} = $value; 
		}
                else {
			my $num = substr($_,1,1);
	                $params_hash{$num} = 'true'; 
                }
	}
	return %params_hash;

};

sub get_parameter_hash{
	my $self = shift;
	my @meta_array = Genome::Model::Tools::Speedseq::Sv-> _tool_param_metas();
	
	my %library = ();
	
	foreach my $meta (@meta_array){
		$library{$meta->tool_param_name} = $meta->property_name; 
	}
	return %library;
}
