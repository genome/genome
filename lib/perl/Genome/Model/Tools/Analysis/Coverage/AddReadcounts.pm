package Genome::Model::Tools::Analysis::Coverage::AddReadcounts;
use strict;
use Genome;
use IO::File;
use warnings;


class Genome::Model::Tools::Analysis::Coverage::AddReadcounts{
    is => 'Command',
    has => [
	bam_files => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'comma-separated list of bam files to grab readcounts from. Output columns will be appended in this order',
	},

	variant_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'File containing snvs in annotation format (1-based, first 5-cols =  [chr, st, sp, ref, var])',
	},

        output_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output file will be indentical to the input file with readcounts appended as the last two columns',
        },

        genome_build => {
            is => 'String',
            is_optional => 0,
	    doc => 'takes either a string describing the genome build (one of 36, 37, mm9, mus37, mus37wOSK) or a path to the genome fasta file',
        },

        min_quality_score => {
            is => 'Integer',
            is_optional => 1,
	    doc => 'minimum mapping quality of reads to be considered',
            default => '1',
        },

        chrom => {
            is => 'String',
            is_optional => 1,
	    doc => 'only process this chromosome.  Useful for enormous files',
        },

        min_depth  => {
            is => 'Integer',
            is_optional => 1,
	    doc => 'minimum depth required for a site to be reported',
        },

        max_depth => {
            is => 'Integer',
            is_optional => 1,
	    doc => 'maximum depth allowed for a site to be reported',
        },

        min_vaf => {
            is => 'Integer',
            is_optional => 1,
	    doc => 'minimum variant allele frequency required for a site to be reported (0-100)',
        },

        max_vaf => { 
            is => 'Integer',
            is_optional => 1,
	    doc => 'maximum variant allele frequency allowed for a site to be reported (0-100)',
        },

        indel_size_limit => {
            is => 'Integer',
            is_optional => 1,
	    doc => 'maximum indel size to grab readcounts for. (The larger the indel, the more skewed the readcounts due ot mapping problems)',
            default => 4,
        },

        header_prefixes => { 
            is => 'String',
            is_optional => 1,
	    doc => 'Comma-separated list - if the file has a header, three column titles get added for each bam ("ref_count","var_count","VAF"). This specifies a prefix for those columns. (i.e.  "Normal" will lead to "Normal_ref_count","Normal_var_count","Normal_VAF").',
        },


        ]
};

sub help_brief {
    "get readcounts. make pretty. append to anno file"
}

sub help_detail {
    "get readcounts. make pretty. append to anno file"
}



sub execute {
    my $self = shift;
    my $bam_files = $self->bam_files;
    my $variant_file = $self->variant_file;
    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $min_quality_score = $self->min_quality_score;

    my $min_vaf = $self->min_vaf;
    my $max_vaf = $self->max_vaf;
    my $min_depth = $self->min_depth;
    my $max_depth = $self->max_depth;
    my $indel_size_limit = $self->indel_size_limit;

    my $chrom = $self->chrom;
    my @header_prefixes;
    my @bams = split(",",$bam_files);
    if(defined($self->header_prefixes)){
        @header_prefixes = split(",",$self->header_prefixes);
    }

    my $fasta;
    if ($genome_build eq "36") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif( ($genome_build eq "37") || ($genome_build eq "37lite") ){
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif ($genome_build eq "mus37") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-mouse-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } elsif ($genome_build eq "mus37wOSK") {
        $fasta = "/gscmnt/sata135/info/medseq/dlarson/iPS_analysis/lentiviral_reference/mousebuild37_plus_lentivirus.fa";
    } elsif ($genome_build eq "mm9") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "UCSC-mouse-buildmm9");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } elsif (-e $genome_build ) {
        $fasta = $genome_build;
    } else {
        die ("invalid genome build or fasta path: $genome_build\n");
    }



    #create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    unless($tempdir) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }


    #set some defaults for undefined parameters
    unless(defined($min_depth)){
        $min_depth = "NA";
    }
    unless(defined($max_depth)){
        $max_depth = "NA";
    }
    unless(defined($min_vaf)){
        $min_vaf = "NA";
    }
    unless(defined($max_vaf)){
        $max_vaf = "NA";
    }
    unless(defined($chrom)){
        $chrom = "all";
    }
   

    my $prefix = 1;
    for my $bam (@bams){
        #run bam-readcount, stick the files in the tempdir
        my $cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
            bam_file => $bam,
            output_file =>  "$tempdir/$prefix.rcfile",
            variant_file => $variant_file,
            genome_build => $genome_build, 
            chrom => $chrom,
            min_depth  => $min_depth,
            max_depth => $max_depth,
            min_vaf => $min_vaf,
            max_vaf => $max_vaf,
            indel_size_limit => $indel_size_limit,
            );
        unless ($cmd->execute) {
            die "Bam-readcount failed";
        }

        my %readcounts;
        #read in the bam-readcount file  and hash each count by position
        my $inFh2 = IO::File->new( "$tempdir/$prefix.rcfile" ) || die "can't open file: $tempdir/$prefix.rcfile\n";
        while( my $line = $inFh2->getline )
        {
            chomp($line);
            my ($chr, $pos, $ref, $var, $refcount, $varcount, $vaf,) = split("\t",$line);
            $ref = uc($ref);
            $var = uc($var);
            
            my $key = "$chr:$pos:$ref:$var";
            
            #for each base at that pos
            $readcounts{$key} = "$refcount:$varcount:$vaf";
        }


        #prep the output file
        open(OUTFILE,">$tempdir/$prefix.output") || die "can't open output for writing\n";

        #read in all the snvs and hash both the ref and var allele by position
        my $inFh;
        if($prefix == 1){
            $inFh = IO::File->new( $variant_file ) || die "can't open file2\n";
        } else {
           $inFh = IO::File->new( "$tempdir/" . ($prefix-1) . ".output" ) || die "can't open file3\n";
        }

        my $count = 0;
        while( my $sline = $inFh->getline )
        {
            chomp($sline);
            
            if($count == 0){ #check for header
                if($sline =~ /^(#|Hugo_Symbol|Chrom|chromosome)/i){
                    #good header match                    
                    if(defined($header_prefixes[$prefix-1])){
                        my $pre = $header_prefixes[$prefix-1];
                        print OUTFILE join("\t",($sline,$pre . "_ref_count", $pre . "_var_count", $pre ."_VAF")) . "\n";
                    } else {
                        print OUTFILE join("\t",($sline,"ref_count","var_count","VAF")) . "\n";
                    }
                    next;
                }
            };
            $count++;
            
            
            my ($chr, $st, $sp, $ref, $var, @rest) = split("\t",$sline);

            $ref =~ s/0/-/g;
            $ref =~ s/\*/-/g;
            $var =~ s/0/-/g;
            $var =~ s/\*/-/g;
            $ref = uc($ref);
            $var = uc($var); 
            
            my $key = $chr . ":" . $st . ":" . $ref . ":" . $var;                
            
            #get the readcount information at this position
            if(exists($readcounts{$key})){
                my ($rcnt,$vcnt,$vaf) = split(/:/,$readcounts{$key});
                $sline = $sline . "\t$rcnt\t$vcnt\t$vaf\n";
            } else {
                print STDERR "didn't find variant $key in bam-readcount output, skipping\n";
            }
            print OUTFILE $sline;
            
        }
        close(OUTFILE);
        $prefix++;
    }

    my $last = $prefix - 1;
    `cp $tempdir/$last.output $output_file`;

    return 1;
}
