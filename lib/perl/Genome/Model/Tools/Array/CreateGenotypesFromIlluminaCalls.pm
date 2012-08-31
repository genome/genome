package Genome::Model::Tools::Array::CreateGenotypesFromIlluminaCalls;

#rlong: This is under development, but seems to be functional

use strict;
use warnings;

use Genome;
use Command;
use Bio::DB::Fasta;
use Text::CSV_XS;
use Sort::Naturally qw( nsort );

class Genome::Model::Tools::Array::CreateGenotypesFromIlluminaCalls {
    is => 'Command',
    has => [
    sample_name =>
    {
        type => 'String',
		is_optional => 0,
		doc => "Name of the sample which is being analyzed",
    },
    call_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => "File containing the calls",
    },
    illumina_manifest_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'File from Illumina indicate genotypes and stranding information',
    },
	output_path =>
	{
		type => 'String',
		is_optional => 0,
		doc => 'Path to the desired location of the output',
	},
	allowed_mismatches =>
	{
		type => 'Number',
		is_optional => 1,
		doc => 'Allowed mismatches, default is 2',
	},
    ]
};


sub execute {
    my $self=shift;
	print "allowed_mismatches is ";

	my $allowed_mismatches;

	if (defined($self->allowed_mismatches)) {
		$allowed_mismatches = $self->allowed_mismatches;
		print "set to ".$allowed_mismatches."\n";
	} else {
		$allowed_mismatches = 2;
		print "set to the default value of 2\n";
	}

    #TODO Some basic file checks
    $| = 1;

    my $call_href = $self->create_call_file_hash;
    my %probes = %$call_href;
    unless($call_href) {
        $self->error_message("Unable to parse genotype file and create call file hash");
        return;
    }
    $self->status_message("Finished creating call file hash");


    #Convert each call into a genotype and write to a new file by chromosome and position
    #create file handles to write out each sample


    #read through the manifest file
    my $bead_fh = new IO::File $self->illumina_manifest_file,"r";
    unless (defined($bead_fh)) {
        $self->error_message("unable to open file ".$self->illumina_manifest_file." for reading.");
	die $self->error_message;
    }

    #get past headers
    while (my $line = $bead_fh->getline) {
        last if $line =~ /IlmnID/;
    }

    #open output file
    my $outfile = $self->output_path."/".$self->sample_name.".genotype";
    my $out_fh = new IO::File $outfile,"w";
    unless (defined($out_fh)) {
        $self->error_message("unable to open file ".$outfile." for reading.");
        die $self->error_message;
    }

    #check which strand probe sequence is on, and write output
    while (my $line = $bead_fh->getline) {
        my ($IlmnID,$probe_name,$IlmnStrand,$SNP,$AddressA_ID,$AlleleA_ProbeSeq,$AddressB_ID,$AlleleB_ProbeSeq,$GenomeBuild,$chr,$pos,$Ploidy,$Species,$Source,$SourceVersion,$SourceStrand,$SourceSeq,$TopGenomicSeq,$BeadSetID,$Intensity_Only,$Exp_Clusters,$CNV_Probe) = split ",",$line;

        #if probe was in dataset
        if (exists $probes{$probe_name}) {

            #Check for a forward strand match
	    (my $source_seq = $SourceSeq) =~ s/(\w+)\[\w\/\w\]\w+/\U$1/;
	    my $probe_length = length($source_seq);
	    if ($chr eq "Mt") { $chr = "MT" };
	    my $ref_seq = get_ref_base($chr,$pos-$probe_length,$pos-1);
	    
	    unless ($ref_seq) {
	        print "Ref seq not found for probe $probe_name, chr $chr, near pos $pos.\n";
	        next;
	    }

	    #if the sequence matches the ref, then probe is on the forward strand (print it out as is)
	    if ($source_seq eq $ref_seq) {
	        $out_fh->print("$chr\t$pos\t$probes{$probe_name}\n");
	        print "Forward-strand match.\n";
	        next;
	    }

	    else {
	        #last attempt - allowing 1 mismatch in the strand vs. the reference
	        my $mismatches;
	        for (my $i = 0; $i < $probe_length; $i++) {
	            my $ref_base = substr($ref_seq,$i,1);
	            my $probe_base = substr($source_seq,$i,1);
	            if ($ref_base eq $probe_base) {next;}
	            else { $mismatches++ };
	        }
	        if ($mismatches <= $allowed_mismatches) {
	            $out_fh->print("$chr\t$pos\t$probes{$probe_name}\n");
	            print "Forward-strand match with $mismatches mismatches.\n";
	            next;
	        }
	    }


	    #No forward strand match. Now, reverse complement the second half of the manifest sequence. If this matches the ref seq, the probe was on rev. strand, so print out rev_comp of called genotypes
	    print "Checking reverse strand for probe $probe_name.\n";
	    ($source_seq = $SourceSeq) =~ s/\w+\[\w\/\w\](\w+)/\U$1/;
	    $probe_length = length($source_seq);
	    if ($probe_name !~ /^rs/) {
	        $pos = $pos - $probe_length*2;
	    }
	    $ref_seq = get_ref_base($chr,$pos-$probe_length,$pos-1);
	    my $rc_source_seq = rev_comp($source_seq);
	    if ($rc_source_seq eq $ref_seq) {
	        $probes{$probe_name} = rev_comp($probes{$probe_name});
	        $out_fh->print("$chr\t$pos\t$probes{$probe_name}\n");
	        print "Reverse-strand match.\n";
	    }
	    else {
	        #last attempt - allowing 1 mismatch in the strand vs. the reference
	        my $mismatches;
	        for (my $i = 0; $i < $probe_length; $i++) {
	            my $ref_base = substr($ref_seq,$i,1);
	            my $probe_base = substr($source_seq,$i,1);
	            if ($ref_base eq $probe_base) {next;}
	            else { $mismatches++ };
	        }
	        if ($mismatches <= $allowed_mismatches) {
	            $probes{$probe_name} = rev_comp($probes{$probe_name});
	            $out_fh->print("$chr\t$pos\t$probes{$probe_name}\n");
	            print "Reverse-strand match with $mismatches mismatches.\n";
	        }
	        else {
	            print "Neither reverse or forward strand match for probe $probe_name, with second half sequence $source_seq. Reverse complement is $rc_source_seq. Meanwhile, ref seq in this case from chr $chr near pos $pos is $ref_seq.\n";
	        }
	    }
	}
	else {
	    print "Probe $probe_name not in dataset.\n";
	}
    }#end, reading through bead manifest file




    $bead_fh->close;
    $out_fh->close;



    return 1;
}

sub create_call_file_hash {
    my $self = shift;
    my $file = $self->call_file;
    #my %call_hash;
	my %probes;
    my $call_fh = new IO::File "$file", "r";
    unless(defined($call_fh)) {
        $self->error_message("Unable to open $file for reading");
        return;
    }

	#get past headers
	while (my $line = $call_fh->getline) {
		last if $line =~ /Allele1/;
	}

    if($call_fh->eof) {
        $self->error_message("Unexpected file format. Encountered end of file before expected");
        return;
    }

	#populate a hash of probe names and results
	while (my $line = $call_fh->getline) {
		my ($sample,$snp,$allele1,$allele2) = split "\t",$line;
		next if $allele1 eq "-";
		next if $allele2 eq "-";
		    $probes{$snp} = "$allele1$allele2";
	}
	$call_fh->close;
	return \%probes;
}

sub get_ref_base {
    my ($chr,$start,$stop) = @_;
    my $ref_dir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    my $refdb = Bio::DB::Fasta->new($ref_dir);
    my $seq = $refdb->seq($chr, $start => $stop);
    $seq =~ s/([\S]+)/\U$1/;
    return $seq;
}

sub rev_comp {
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/actgACTG/TGACTGAC/;
    return $seq;
}


1;

sub help_brief {
    "Converts an Illumina Beadstudio Forward Strand Report file to genotype file";
}

