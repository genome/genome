package Genome::Model::Tools::CopyNumber::CreateIlluminaCnFiles;

use warnings;
use strict;
use Genome;
use Cwd;
use IO::File;
use Getopt::Long;
use Data::Dumper;

class Genome::Model::Tools::CopyNumber::CreateIlluminaCnFiles {
    is => 'Command',
    has => [
    output_dir => {
        is => 'String',
        is_optional => 1,
        default => getcwd(),
        doc => 'Directory containing output files and directories.',
    },
    annotation_file => {
        is => 'String',
        is_optional => 0,
        doc => 'Illumina annotation.txt file',
    },
    copy_number_file => {
        is => 'String',
        is_optional => 0,
        doc => 'Illumina copy number output file, typically "...pairedcopynumber.txt".',
    },
    ]
};

sub help_brief {
    "generate copynumber data file for each sample and map.csv file for grouped data"
}

sub help_detail {
    "Use this script to process Illumina SNP array data. It will generate a map.csv file mapping out positions of data files for all samples of a grouped set of data, and it will generate a folder /CN/ in the output dir containing a copynumber data file for each sample."
}
    
sub execute {
    my $self = shift;

    #process input arguments
    my $output_dir = $self->output_dir;
    `mkdir $output_dir` unless (-e "$output_dir");
    my $anno_file = $self->annotation_file;
    my $anno_fh = new IO::File $anno_file,"r";
    my $cn_file = $self->copy_number_file;
    
    #read annotation file
    my %pos_hash;
    while (my $line = $anno_fh->getline) {
        next if ($line =~ /^Name/);
        #chomp $line;
        my ($snp_id,$index,$chr,$pos) = split /\t/, $line;

        next if $chr eq "XY";
        $pos_hash{$snp_id}{chr} = $chr;
        $pos_hash{$snp_id}{pos} = $pos;
    }

    #copy number file is usually in dos format - fix this!
    `dos2unix -bd $cn_file`;

    #open copy number file and parse headers
    my $cn_fh = new IO::File $cn_file,"r";
    my $header = $cn_fh->getline;
    chomp($header);
    #Sample ID       Gender     Sample Name     Sample Group    Call Rate       Phenotype
    my ($Sample_id,$Gender,$Sample_name,$Sample_group,$Call_rate,$phenotype,@snp_ids) = split /\t/, $header;

    #write unsorted SNP file
    my $unsorted_SNP_fh = new IO::File "$output_dir/SNP.csv.unsorted","w";
    my %map_hash;
    for my $id (@snp_ids) {
        if (defined $pos_hash{$id}{chr} && defined $pos_hash{$id}{pos}) {
            my $chr=$pos_hash{$id}{chr};
            my $pos=$pos_hash{$id}{pos};
            print $unsorted_SNP_fh "$chr\t$pos\n";
            my $ch=$chr;
            $ch=23 if ($ch eq "MT");
            $ch=24 if ($ch eq "X");
            $ch=25 if ($ch eq "Y");
            $map_hash{$ch}{$pos}=$id;
        }
    }
    $unsorted_SNP_fh->close;

    #write map.csv
    my $map_fh = new IO::File "$output_dir/map.csv","w";
    print $map_fh "CHR\tPOS\n";
    for my $ch (sort numer keys %map_hash) {
        for my $pos (sort numer keys %{$map_hash{$ch}}) {
            my $CH=$ch;
            $CH="MT" if ($CH eq "23");
            $CH="X" if ($CH eq "24");
            $CH="Y" if ($CH eq "25");
            print $map_fh "$CH\t$pos\n";
        }
    }
    $map_fh->close;


    my %cn_hash;
    `mkdir $output_dir/CN` unless (-e "$output_dir/CN");
    while (my $line = $cn_fh->getline) { #rest of cn lines will generate log2 CN data (sorted)
        chomp $line;
        my ($sample_id,$Gender,$Sample_name,$Sample_group,$Call_rate,$phenotype,@cns) = split /\t/, $line;
        #my ($sample_id,$sample_pair,$gender,$sample_group,@cns) = split /\t/, $line;

        if ($sample_id=~/\s/)
        {
            my @parts = split /\s+/, $sample_id;
            $sample_id=join "_",@parts;
        }
        print STDOUT "Working on sample $sample_id.\n";
        my $cn_file="$output_dir/CN/$sample_id.log2";
        my $cn_fh = new IO::File $cn_file,"w";
        print $cn_fh "Log2_Copy_Number_Data_For_".$sample_id."\n";
        for(my $i=0; $i<=$#cns; $i++)
        {
            my $snp_id=$snp_ids[$i];
            my $chr=$pos_hash{$snp_id}{chr};
            my $pos=$pos_hash{$snp_id}{pos};
            if ($chr && $pos && $cns[$i])
            {
                $cn_hash{$chr}{$pos}=$cns[$i];
            }
            else
            {
                print "please check - snp_id: $snp_id, chr: $chr, pos: $pos, cn: $cns[$i]\n";
            }
        }

        for my $ch (1..22,"MT","X","Y")
        {
            for my $pos (sort numer keys %{$cn_hash{$ch}})
            {
                my $cn=$cn_hash{$ch}{$pos};
                my $log_cn;
                #print "$ch\t$pos\t$cn\n";
                if ($cn=~/^\d+/)
                {
                    if ($cn==0)
                    {
                        $log_cn=-10;
                    }
                    else
                    {
                        $log_cn=log($cn/2)/log(2); #copy_number=2*2**(log2intensity)
                    }
                    printf $cn_fh "%-8.4f\n", $log_cn;
                }
                else
                {
                    $log_cn="NA";
                    print $cn_fh "$log_cn\n";
                }
            }

        }

    }

return 1;
}

sub numer { 
    $a <=> $b
}

1;
