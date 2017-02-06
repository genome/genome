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
            doc => 'takes either a string describing the genome build (one of 36, 37, mm9, mus37) or a path to the genome fasta file',
        },
        min_quality_score => {
            is => 'Integer',
            is_optional => 1,
            doc => 'minimum mapping quality of reads to be considered',
            default => '1',
        },
        min_base_quality => {
            is => 'Integer',
            is_optional => 1,
            doc => 'minimum base quality of bases in reads to be considered',
            default => '0',
        },
        chrom => {
            is => 'String',
            is_optional => 1,
            doc => 'only process this chromosome.  Useful for enormous files',
        },
        min_depth  => {
            is => 'Integer',
            is_optional => 1,
            doc => 'minimum depth required for a site to be reported. Only use with a single BAM.',
        },
        max_depth => {
            is => 'Integer',
            is_optional => 1,
            doc => 'maximum depth allowed for a site to be reported. Only use with a single BAM.',
        },
        min_vaf => {
            is => 'Integer',
            is_optional => 1,
            doc => 'minimum variant allele frequency required for a site to be reported (0-100). Only use with a single BAM.',
        },
        max_vaf => {
            is => 'Integer',
            is_optional => 1,
            doc => 'maximum variant allele frequency allowed for a site to be reported (0-100). Only use with a single BAM.',,
        },
        indel_size_limit => {
            is => 'Integer',
            is_optional => 1,
            doc => 'maximum indel size to grab readcounts for. (The larger the indel, the more skewed the readcounts due to mapping problems)',
            default => 4,
        },
        header_prefixes => {
            is => 'String',
            is_optional => 1,
            doc => 'Comma-separated list - if the file has a header, three column titles get added for each bam ("ref_count","var_count","VAF"). This specifies a prefix for those columns. (i.e.  "Normal" will lead to "Normal_ref_count","Normal_var_count","Normal_VAF").',
        },
        per_library => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'whether or not to report counts on a per-library basis',
            default => 0,
        },
        bam_readcount_version => {
            is => 'String',
            doc => 'version of bam-readcount to use',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    "get readcounts. make pretty. append to anno file"
}

sub help_detail {
    "This expects input in old TGI annotator format. It appends reference, variant counts and VAF (as a percentage) to the end of each input line."
}



sub execute {
    my $self = shift;
    my @bams = split(",", $self->bam_files);
    my $variant_file = $self->variant_file;
    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $min_quality_score = $self->min_quality_score;
    my $min_base_quality = $self->min_base_quality;

    # These options result in sites failing the filters to be removed
    # from the output of Genome::Model::Tools::Analysis::Coverage::BamReadcount
    # This is appropriate for single samples, but not for multiple. We will check
    # and die with an appropriate error if someone tries to do the wrong thing

    my $min_vaf = $self->min_vaf;
    my $max_vaf = $self->max_vaf;
    my $min_depth = $self->min_depth;
    my $max_depth = $self->max_depth;

    if(scalar(@bams) > 1 && grep { defined $_ } ($min_vaf, $max_vaf, $min_depth, $max_depth)) {
        die "Cannot filter on VAF and/or Depth if gathering results from multiple BAM files";
    }

    my $indel_size_limit = $self->indel_size_limit;

    my $chrom = $self->chrom;
    my @header_prefixes;
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
    } elsif ($genome_build eq "mm9") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(id => "107494762");
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
        my %params = (
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
            min_mapping_quality => $min_quality_score,
            min_base_quality => $min_base_quality,
            per_library => $self->per_library,
            );
        if ($self->bam_readcount_version) {
            $params{bam_readcount_version} = $self->bam_readcount_version;
        }
        my $cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(%params);
        unless ($cmd->execute) {
            die "Bam-readcount failed";
        }

        my %readcounts;
        my @per_lib_headers;
        #read in the bam-readcount file  and hash each count by position
        my $inFh2 = IO::File->new( "$tempdir/$prefix.rcfile" ) || die "can't open file: $tempdir/$prefix.rcfile\n";
        while( my $line = $inFh2->getline )
        {
            chomp($line);
            if($line =~ /^#chr/) {
                $line =~ s/^#//;
                @per_lib_headers = split("\t",$line);
                next;
            }
            my ($chr, $pos, $ref, $var, @counts,) = split("\t",$line);
            $ref = uc($ref);
            $var = uc($var);

            my $key = "$chr:$pos:$ref:$var";

            #for each base at that pos
            $readcounts{$key} = join(":", @counts);
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
        my @count_headers = qw( ref_count var_count VAF);
        if($self->per_library) {
            @count_headers = @per_lib_headers[4..$#per_lib_headers];
        }
        while( my $sline = $inFh->getline )
        {
            chomp($sline);

            if($count == 0){ #check for header
                if($sline =~ /^(#|Hugo_Symbol|Chrom|chromosome|chr\s)/i) {
                    #good header match
                    if(defined($header_prefixes[$prefix-1])){
                        my $pre = $header_prefixes[$prefix-1];
                        print OUTFILE join("\t",($sline, map { $pre . "_" . $_ } @count_headers)) . "\n";
                    } else {
                        print OUTFILE join("\t",($sline,@count_headers)) . "\n";
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
                my @counts = split(/:/,$readcounts{$key});
                $sline = $sline . "\t" . join("\t", @counts);
            } else {
                print STDERR "didn't find variant $key in bam-readcount output, skipping\n";
            }
            print OUTFILE $sline, "\n";

        }
        close(OUTFILE);
        $prefix++;
    }

    my $last = $prefix - 1;
    my $cmd = "cp $tempdir/$last.output '" . $output_file . "'";
    `$cmd`;

    return 1;
}
