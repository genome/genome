package Genome::Model::Tools::Annotate::AddRsid;

use strict;
use warnings;
use Genome;


class Genome::Model::Tools::Annotate::AddRsid {
    is => ['Command'],
    has_input => [
        anno_file => {
            doc => 'The name of the annotation file to use. default_value=',
        },
        vcf_file => {
            doc => 'The vcf file containing RSid information',
        },
	output_file => {
	    
	},
       
    ],
};

sub execute {

    $DB::single=1;
    my $self = shift;

    my $anno_file = $self->anno_file;
    my $vcf_file = $self->vcf_file;
    my $output_file = $self->output_file;
    die "$output_file already exists.  Aborting...\n" if(-e $output_file);

    my $RSid = store_RSid($vcf_file);

    my $out_fh = Genome::Sys->open_file_for_writing($output_file);
    my $anno_fh = Genome::Sys->open_file_for_reading($anno_file);
    while(<$anno_fh>) {
	next if(/^chromosome/);
	chomp;
	my @list = split(/\t/,$_);
	my $k = join("_",@list[0,1]);
	if(exists($RSid->{$k}->{$list[4]})) {
	    push(@list,$RSid->{$k}->{$list[4]});
	}else {
	    push(@list,"-");
	}
	my $str = join("\t",@list);
	$out_fh->print("$str\n");

    }
    $anno_fh->close;
    $out_fh->close;

    
    return 1;
}

sub store_RSid {

    my $vcf = shift;

    my $vcf_fh = Genome::Sys->open_gzip_file_for_reading($vcf);
    return unless(defined($vcf_fh)); 

    my $RSid={};
    while(<$vcf_fh>) {
	next if(/^\#/);
	chomp;
	my ($chr,$pos,$rsID,$ref,$var,@rest) = split(/\t/,$_);
	next if($rsID eq '.'); #skip if RSid not defined

	my $key = join("_",($chr,$pos));
	my @var_alleles = split(/,/,$var);
	my $RSid_var_allele;
	if(@var_alleles > 1) { #multiple variant alleles
	    my $INFO=$rest[2];
	    my ($dbSNPinfo) = $INFO =~ /dbSNPBuildID=([0-9,.]+)/;
	    if($dbSNPinfo !~ /\./) {
		#print STDERR "Error, multiple dbSNPbuildVersion found for $chr,$pos,$ref:$var,$rsID,$rest[2]\n";
	    }
	    my @dbSNPids = split(/,/,$dbSNPinfo);
	    for(my $i=0;$i<@dbSNPids;$i++) {
		if($dbSNPids[$i] =~ /^\d+$/) {
		    $RSid_var_allele = $var_alleles[$i];
		    $RSid->{$key}->{$RSid_var_allele}=$rsID;
		}
	    }
	}else { #just 1 variant allele
	    $RSid_var_allele = $var;
	    $RSid->{$key}->{$RSid_var_allele}=$rsID;
	}
	

    }
    $vcf_fh->close();

    print STDERR "done parsing VCF file\n";
    return $RSid;

}

 
1;
