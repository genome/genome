#!/usr/bin/env genome-perl
#convert from annot file format to VCF format
package Genome::Model::Tools::Tcga::ConvertMfToVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use Time::localtime;

class Genome::Model::Tools::Tcga::ConvertMfToVcf {
	is => 'Genome::Model::Tools::Tcga',
	has => [
		mf_file => {
			type => 'String',
			doc => 'input MF file, default in annotation format',
		},
		output_file => {
			type => 'String',
			doc => 'output file in VCF file format',
		},
	],
	has_optional => [
                _not_overwrite => {
                        type => 'BOOLEAN',
                        default => 0,
                },
		_mf_fh => {
			type => 'SCALAR',
		},
		_out_fh => {
			type => 'SCALAR',
		}
	],
};

sub execute {
	my $self = shift;
	my $mf_file = $self->mf_file;
	my $mf_fh = Genome::Sys->open_file_for_reading($mf_file) or return;
	my $datetime_string = ctime(stat($mf_file)->mtime); # read the file time
	$self->_mf_fh($mf_fh);
	my @VCF;
	while(my $line = $mf_fh->getline) {
		if($line =~ /^#/){
			next;
		}
		chomp $line;
		my $vcf;
		my @line_ = split /\s+/, $line;
		$vcf->{line} = $line;
		$vcf->{chr} = $line_[0];
		$vcf->{pos} = $line_[1];
		$vcf->{ID} = ".";
		$vcf->{REF} = $line_[3];
		$vcf->{ALT} = $line_[4];
                $vcf->{QUAL} = ".";
                $vcf->{QUAL} = $line_[$#line_-3] if($line_[$#line_-3] =~ /^\d+$/);
		$vcf->{FILTER} = "PASS";
		my $size = abs($line_[2] - $line_[1] + 1);
		$vcf->{INFO} = "SZ=" . $size . ";TP=" . $line_[5] . ";GT=" . $line_[$#line_-2] . ";TR=". $line_[$#line_-1]. ";CD=" . $line_[$#line_];
		$vcf->{REF} = "." if($vcf->{REF} eq "0");
		$vcf->{ALT} = "." if($vcf->{ALT} eq "0");
		push @VCF, $vcf;
	}

	my $out_file = $self->output_file;
        `rm -f $out_file` if(! $self->_not_overwrite);
	my $out_fh = Genome::Sys->open_file_for_writing($out_file) or return;
	$self->_out_fh($out_fh);

	$out_fh->print("##fileformat=VCFv4.0\n");
	$out_fh->print("##fileDate=$datetime_string\n");
	$out_fh->print("##source=GSCprogram\n");
	$out_fh->print("##reference=NCBI36\n");
	$out_fh->print("##phasing=partial\n");
	$out_fh->print("##INFO=<ID=SZ,Number=1,Type=Integer,Description=\"Size\">\n");
	$out_fh->print("##INFO=<ID=TP,Number=1,Type=String,Description=\"Type\">\n");
        $out_fh->print("##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        $out_fh->print("##INFO=<ID=TR,Number=1,Type=Interger,Description=\"Tier from 1 to 4\">\n");
        $out_fh->print("##INFO=<ID=CD,Number=1,Type=String,Description=\"Confidence of High or Low\">\n");
	# no FILTER, no FORMAT
	$out_fh->print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

	for my $vcf (@VCF) {
		$out_fh->print("$vcf->{chr}\t$vcf->{pos}\t$vcf->{ID}\t$vcf->{REF}\t$vcf->{ALT}\t$vcf->{QUAL}\t$vcf->{FILTER}\t$vcf->{INFO}\n");
	}
	$self->_mf_fh->close;
	$self->_out_fh->close;
	return 1;
}

1;


