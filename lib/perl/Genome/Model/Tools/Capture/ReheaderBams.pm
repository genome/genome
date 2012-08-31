
package Genome::Model::Tools::Capture::ReheaderBams;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReheaderBams - Input a bam, output a bam with new header based on options
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	02/04/2011 by W.S.
#	MODIFIED:	02/04/2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Capture::ReheaderBams {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		input_bam_file	=> { is => 'Text', doc => "Bam File" , is_optional => 0},
		output_bam_file	=> { is => 'Text', doc => "Bam File" , is_optional => 0},
		old_as	=> { is => 'Text', doc => "Original AS Field to Change, ex: GRCh37-lite-build37" , is_optional => 1},
		old_ur	=> { is => 'Text', doc => "Original UR Field to Change, ex: ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz" , is_optional => 1},
#		old_sn	=> { is => 'Text', doc => "Original SN Field to Change" , is_optional => 1},
#		old_ln	=> { is => 'Text', doc => "Original LN Field to Change" , is_optional => 1},
#		old_m5	=> { is => 'Text', doc => "Original M5 Field to Change" , is_optional => 1},
#		old_sp	=> { is => 'Text', doc => "Original SP Field to Change" , is_optional => 1},
		new_as	=> { is => 'Text', doc => "Changes to New AS Field" , is_optional => 1},
		new_ur	=> { is => 'Text', doc => "Changes to New UR Field" , is_optional => 1},
#		new_sn	=> { is => 'Text', doc => "Changes to New SN Field" , is_optional => 1},
#		new_ln	=> { is => 'Text', doc => "Changes to New LN Field" , is_optional => 1},
#		new_m5	=> { is => 'Text', doc => "Changes to New M5 Field" , is_optional => 1},
#		new_sp	=> { is => 'Text', doc => "Changes to New SP Field" , is_optional => 1},

	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Input a bam, output a bam with new header based on options"                 
}

sub help_synopsis {
    return <<EOS
Tool to reheader the bam files, used primarily for TCGA bams
EXAMPLE:	gmt capture reheader-bams
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $old_AS = "AS:" . $self->old_as if $self->old_as;
	my $old_UR = "UR:" . $self->old_ur if $self->old_ur;
#	my $old_SN = "SN:" . $self->old_sn if $self->old_sn;
#	my $old_LN = "LN:" . $self->old_ln if $self->old_ln;
#	my $old_M5 = "M5:" . $self->old_m5 if $self->old_m5;
#	my $old_SP = "SP:" . $self->old_sp if $self->old_sp;
	my $new_AS = "AS:" . $self->new_as if $self->new_as;
	my $new_UR = "UR:" . $self->new_ur if $self->new_ur;
#	my $new_SN = "SN:" . $self->new_sn if $self->new_sn;
#	my $new_LN = "LN:" . $self->new_ln if $self->new_ln;
#	my $new_M5 = "M5:" . $self->new_m5 if $self->new_m5;
#	my $new_SP = "SP:" . $self->new_sp if $self->new_sp;

	if ($self->old_as || $self->new_as) {
		unless ($self->old_as && $self->new_as) {
			die "Cannot have old_AS $old_AS or new_AS $new_AS, must have both or neither defined";
		}
	}
	if ($self->old_ur || $self->new_ur) {
		unless ($self->old_ur && $self->new_ur) {
			die "Cannot have old_ur $old_UR or new_ur $new_UR, must have both or neither defined";
		}
	}
	
	my $input_bam_file = $self->input_bam_file;
	my $output_bam_file = $self->output_bam_file;

	my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
	unless($tfh) {
		$self->error_message("Unable to create temporary file $!");
		die;
	}

	my ($tfh2,$temp_path2) = Genome::Sys->create_temp_file;
	unless($tfh2) {
		$self->error_message("Unable to create temporary file $!");
		die;
	}
	open(OUTFILE1, ">$temp_path2") or die "Can't open outfile: $!\n";

	my $cmd1 = "samtools view -H $input_bam_file > $temp_path";

	my $return = Genome::Sys->shellcmd(
		cmd => "$cmd1",
		output_files => [$temp_path],
		skip_if_output_is_present => 0,
	);
	unless($return) { 
		$self->error_message("Failed to execute samtools header file creation: $cmd1 Returned $return");
		die $self->error_message;
	}

	my $input = new FileHandle ($temp_path);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;

		my $new_line = "";

		for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
		{
			$new_line .= "\t" if($new_line);
			if ($self->old_as && $self->new_as) {
				if($lineContents[$colCounter] && $lineContents[$colCounter] =~ $old_AS)	{
					$lineContents[$colCounter] = $new_AS;
				}
				elsif($lineContents[$colCounter] && $lineContents[$colCounter] =~ 'AS')	{
					print "Warning: $lineContents[$colCounter] does not match $old_AS\n";
				}
			}

			if ($self->old_ur && $self->new_ur) {
				if($lineContents[$colCounter] && $lineContents[$colCounter] eq $old_UR)	{
					$lineContents[$colCounter] = $new_UR;
				}
				elsif($lineContents[$colCounter] && $lineContents[$colCounter] =~ 'UR')	{
					print "Warning: Old path does not match!\n$lineContents[$colCounter] (header) != \n$old_UR (expected)\n";
				}
			}
			
			$new_line .= $lineContents[$colCounter];
		}

		print OUTFILE1 "$new_line\n";
	}
	
	close($input);	
	close(OUTFILE1);

	my $cmd2 = "samtools reheader $temp_path2 $input_bam_file > $output_bam_file";
	$return = Genome::Sys->shellcmd(
		cmd => "$cmd2",
		output_files => [$output_bam_file],
		skip_if_output_is_present => 0,
	);
	unless($return) { 
		$self->error_message("Failed to execute samtools reheader file creation: $cmd2 Returned $return");
		die $self->error_message;
	}
}

1;
