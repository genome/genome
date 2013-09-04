package Genome::Model::Tools::Analysis::Coverage::FormatReadcounts;

use strict;
use Genome;
use IO::File;
use warnings;


class Genome::Model::Tools::Analysis::Coverage::FormatReadcounts{
    is => 'Command',
    has => [
	snv_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'File containing snvs in 1-based, 5-col format (chr, st, sp, var, ref)',
	},

        output_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output file suitable for input into clonality plot',
        },

        readcount_file => {
            is => 'String',
            is_optional => 0,
	    doc => 'file produced by bam-readcount0.4 (or earlier)',
            default => '1',
        },

        ]
};

sub help_brief {
    "given a set of snvs and raw readcounts from bam-readcount0.4, parse out the variant and ref alleles, and calculate the variant allele frequencies"
}

sub help_detail {
    "given a set of snvs and raw readcounts from bam-readcount0.4, parse out the variant and ref alleles, and calculate the variant allele frequencies"
}



sub execute {
    my $self = shift;
    my $snv_file = $self->snv_file;
    my $output_file = $self->output_file;
    my $readcount_file = $self->readcount_file;


    my %refHash;
    my %varHash;

    #read in all the snvs and hash both the ref and var allele by position
    my $inFh = IO::File->new( $snv_file ) || die "can't open file\n";
    while( my $sline = $inFh->getline )
    {
        chomp($sline);
        my @fields = split("\t",$sline);
        $refHash{$fields[0] . "|" . $fields[1]} = $fields[3];
        $varHash{$fields[0] . "|" . $fields[1]} = $fields[4]
    }
    

    #convert iub bases to lists
    sub convertIub{
	my ($base) = @_;
	#deal with cases like "A/T" or "C/T"
	if ($base =~/\//){
	    my @bases=split(/\//,$base);
	    my %baseHash;
	    foreach my $b (@bases){
		my $res = convertIub($b);
		my @bases2 = split(",",$res);
		foreach my $b2 (@bases2){
		    $baseHash{$b2} = 0;
		}
	    }
	    return join(",",keys(%baseHash));
	}

    return join(',', Genome::Info::IUB->iub_to_bases($base)); #TODO: this is completely stupid.  make this not
    };

    sub matchIub{
        my ($allele,$ref,$var) = @_;
        my @variubs = split(",",convertIub(uc($var)));
        my @refiubs = split(",",convertIub(uc($ref)));
        foreach my $i (@variubs){
            unless (grep {$_ eq $i} @refiubs) {
                if ($allele eq $i){
                    return 1;
                }
            }
        }
        return 0;
    }



    #prep the output file
    open(OUTFILE,">$output_file") || die "can't open $output_file for writing\n";
   

    #read in the bam-readcount file
    my $inFh2 = IO::File->new( "$readcount_file" ) || die "can't open file\n";
    while( my $line = $inFh2->getline )
    {
        chomp($line);
        my ($chr, $pos, $ref, $depth, @counts) = split("\t",$line);

        my $ref_count = 0;
        my $var_count = 0;
        my $knownRef;
        my $knownVar;
        my $var_freq = 0;
        
        # skip if it's not in our list of snvs
        next unless (exists($refHash{$chr . "|" . $pos}) && exists($varHash{$chr . "|" . $pos}));

        #for each base at that pos
        foreach my $count_stats (@counts) {
            my ($allele, $count, $mq, $bq) = split /:/, $count_stats;
            
            #look up the snv calls at this position
            $knownRef = $refHash{$chr . "|" . $pos};
            $knownVar = $varHash{$chr . "|" . $pos};

            #handle snvs first
            if($knownRef ne "-" && $knownVar ne "-"){
                # assume that the ref call is ACTG, not iub 
                # (assumption looks valid in my files)
                if ($allele eq $knownRef){
                    $ref_count += $count;
                }
                
                # if this base is included in the IUB code for
                # for the variant, (but doesn't match the ref)
                if (matchIub($allele,$knownRef,$knownVar)){
                    $var_count += $count;
                }
                
                if ($depth ne '0') {
                    $var_freq = $var_count/$depth * 100;
                }            

            } else { #is an indel, skip it
                $ref_count = "NA";
                $var_count = "NA";
                $var_freq = "NA";
            }
        }

        print OUTFILE "$chr\t$pos\t$knownRef\t$knownVar\t$ref_count\t$var_count\t";
        if ($var_freq eq "NA"){
            print OUTFILE $var_freq;
        } else {
            print OUTFILE sprintf("%.2f",$var_freq);
        }
        print OUTFILE "\n";
    }
    close(OUTFILE)
}
