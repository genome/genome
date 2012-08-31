package Genome::Model::Tools::Vcf::VcfSplitSamples;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;

class Genome::Model::Tools::Vcf::VcfSplitSamples {
    is => 'Command',
    has => [
	output_dir => {
	    is => 'Text',
	    is_optional => 0,
	    doc => "Outputs vcf file for each sample into directory",
	},
        vcf_input => {
            is => 'Text',
            is_optional => 0,
            doc => "VCF file containing mutations from multiple samples",
        },
	supplement_filename => {
            is => 'Text',
            is_optional => 1,
            doc => "default is samplename.vcf if you add this option it will be samplename.supplement.vcf",
        },
	],
};


sub help_brief {                            # keep this to just a few words <---
    "Split single multi-sample VCF into individual single-sample VCFs"
}


sub help_synopsis {
<<'HELP';
Split single multi-sample VCF into individual single-sample VCFs
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
<<'HELP';
Split single multi-sample VCF into individual single-sample VCFs
HELP
}

###############

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $vcf_input = $self->vcf_input;
	my $output_dir = $self->output_dir;

	my @headerlines;
	my @samples;
	my %sample_hash;
    my $inFh;
    if(Genome::Sys->_file_type($vcf_input) eq 'gzip') {
        $inFh = Genome::Sys->open_gzip_file_for_reading($vcf_input);
    }
    else {
	    $inFh = IO::File->new($vcf_input) || die "can't open file\n";
    }
	while(my $line = $inFh->getline ) {
        chomp($line);
        if ($line =~ /^##/){
			push(@headerlines,$line);
		}
    	elsif ($line =~ /^#/){
    		if ($line =~ /#CHROM/){
#				$line =~ s/#//;
				my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample_data) = split(/\t/, $line);
				@samples = @sample_data;
				print "First sample is: ".$samples[0]."\n";
				foreach my $sample (@sample_data) {
					my $data = "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample\n";
					$sample_hash{$sample}{$data}++;
#					my $output_file = $output_dir . "/$sample.vcf";
#					my $fh;
#					print $fh$sample "$header\n";
#					print $fh$sample "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample\n";
				}
			}
			else {
				die "Header is not what we expect sample header line to look like. Please reformat your file or modify script to handle your file.\n";
			}
    	}
		else {
			my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample_data) = split(/\t/, $line);
			my $length = @samples;
			my $i = 0;
			while ($i < $length) { #want (length - 1) since it's 0 based.
				my $current_sample = $samples[$i];
				my $current_geno = $sample_data[$i];
				my $data = "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$current_geno\n";
				$sample_hash{$current_sample}{$data}++;
#				print $fh$sample "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t";
				$i++;
			}
		}
	}
	foreach my $sample (sort keys %sample_hash) {
		my $output_file;
		if ($self->supplement_filename) {
			my $supplement_filename = $self->supplement_filename;
			$output_file = $output_dir . "/$sample.$supplement_filename.vcf";
		}
		else {
			$output_file = $output_dir . "/$sample.vcf";
		}
		open(OUTFILE, ">" . $output_file) or die "Can't open outfile: $!\n";
		my $header = join("\n",@headerlines);
		print OUTFILE "$header\n";
		foreach my $dataline (sort keys %{$sample_hash{$sample}}) {
			print OUTFILE "$dataline";
		}
		close (OUTFILE);
	}

	return 1;

}


