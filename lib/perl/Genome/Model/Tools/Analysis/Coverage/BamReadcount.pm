package Genome::Model::Tools::Analysis::Coverage::BamReadcount;
use strict;
use Genome;
use IO::File;
use warnings;
use Genome::Model::Tools::Vcf::Helpers qw/convertIub/;
use Genome::File::BamReadcount::Reader;



class Genome::Model::Tools::Analysis::Coverage::BamReadcount{
    is => 'Command',
    has => [
    bam_file => {
        is => 'String',
        is_optional => 0,
        doc => 'path to the bam file to get readcounts from',
    },

    variant_file => {
        is => 'String',
        is_optional => 0,
        doc => 'File containing snvs in 1-based, 5-col format (chr, st, sp, ref, var)',
    },

    output_file => {
        is => 'String',
        is_optional => 0,
        doc => 'output file (chr, position, ref_count, var_count, var_freq)',
    },

    genome_build => {
        is => 'String',
        is_optional => 0,
        doc => 'takes either a string describing the genome build (one of 36, 37, mm9, mus37, mus37wOSK) or a path to the genome fasta file',
    },

    min_mapping_quality => {
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
        is => 'String',
        is_optional => 1,
        doc => 'minimum depth required for a site to be reported',
    },

    max_depth => {
        is => 'String',
        is_optional => 1,
        doc => 'maximum depth allowed for a site to be reported',
    },

    min_vaf => {
        is => 'String',
        is_optional => 1,
        doc => 'minimum variant allele frequency required for a site to be reported (0-100)',
    },

    max_vaf => {
        is => 'String',
        is_optional => 1,
        doc => 'maximum variant allele frequency allowed for a site to be reported (0-100)',
    },

    indel_size_limit => {
        is => 'Integer',
        is_optional => 1,
        doc => 'maximum indel size to grab readcounts for. (The larger the indel, the more skewed the readcounts due to mapping problems)',
        default => 4,
    },

    count_non_reference_reads => {
        is => 'Boolean',
        is_optional => 1,
        doc => 'if this flag is set, the tool will return the count and frequency of all non-reference reads, not just the frequency of the variant listed. Currently only works on SNVs, will skip indels'
    },

    per_library => {
        is => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => 'generate a file containing per-library counts. Will contain a header.',
    },

    ]
};

sub help_brief {
    "get readcounts. make pretty. output ref, var, vaf"
}

sub help_detail {
    "get readcounts. make pretty"
}


## This process could be done more efficiently (less hashes, etc), but this works for now.

sub execute {
    my $self = shift;
    my $bam_file = $self->bam_file;
    my $variant_file = $self->variant_file;
    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $min_base_quality = $self->min_base_quality;
    my $min_mapping_quality = $self->min_mapping_quality;
    my $indel_size_limit = $self->indel_size_limit;
    my $count_non_reference_reads = $self->count_non_reference_reads;
    if($count_non_reference_reads){
        $indel_size_limit = 0;
    }

    my $min_vaf = $self->min_vaf;
    my $max_vaf = $self->max_vaf;
    my $min_depth = $self->min_depth;
    my $max_depth = $self->max_depth;
    my $chrom = $self->chrom;

    #grab the appropriate fasta file
    my $fasta;
    if ($genome_build eq "36") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif (($genome_build eq "37") || ($genome_build eq "37lite")) {
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



    #--------------------------------------------------
    #convert iub bases to lists

    sub matchIub{
        my ($allele,$ref,$var) = @_;
        my @variubs = split(",",convertIub($var));
        my @refiubs = split(",",convertIub($ref));
        foreach my $i (@variubs){
            unless (grep {$_ eq $i} @refiubs) {
                if ($allele eq $i){
                    return 1;
                }
            }
        }
        return 0;
    }

    sub shouldFilter {
        my ($ref_count, $var_count, $var_freq, $min_depth, $max_depth, $min_vaf, $max_vaf) = @_;
        if($var_freq eq "NA") {
            return 0;
        }
        else {
            if(defined($min_depth) && $min_depth ne "NA") {
                return 1 if(($ref_count + $var_count) < $min_depth);
            }
            if(defined($max_depth) && $max_depth ne "NA") {
                return 1 if(($ref_count + $var_count) > $max_depth);
            }
            if(defined($min_vaf) && $min_vaf ne "NA") {
                return 1 if( $var_freq < $min_vaf );
            }
            if(defined($max_vaf) && $max_vaf ne "NA") {
                return 1 if( $var_freq > $max_vaf)
            }
        }
        return 0;
    }

    sub filterAndPrint{
        my ($chr, $pos, $knownRef, $knownVar, $ref_count, $var_count, $var_freq,
            $min_depth,$max_depth,$min_vaf,$max_vaf,$OUTFILE) = @_;
        #handle the special case where this value is NA, which means don't filter, but
        #pass everything through. This lets AddReadcounts.pm work correctly.
        unless(shouldFilter($ref_count, $var_count, $var_freq, $min_depth, $max_depth, $min_vaf, $max_vaf)) {
            print $OUTFILE "$chr\t$pos\t$knownRef\t$knownVar\t$ref_count\t$var_count\t";
            if ($var_freq eq "NA"){
                print $OUTFILE $var_freq;
            } else {
                print $OUTFILE sprintf("%.2f",$var_freq);
            }
            print $OUTFILE "\n";
        }
    }

    sub printLibs {
        my ($OUTFILE, $chr, $pos, $ref, $var, @counts) = @_;
        print $OUTFILE "$chr\t$pos\t$ref\t$var";
        while(@counts) {
            my ($ref_count, $var_count, $var_freq) = splice @counts, 0, 3;
            print $OUTFILE "\t",join("\t", $ref_count, $var_count),"\t";
            if ($var_freq eq "NA"){
                print $OUTFILE $var_freq;
            } else {
                print $OUTFILE sprintf("%.2f",$var_freq);
            }
        }
        print $OUTFILE "\n";
    }



    #---------------------------

    #create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    unless($tempdir) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }

    #split out the chromosome we're working on, if necessary
    if (defined($chrom) && ($chrom ne "all")){
        my $cmd = "grep \"^" . $chrom . "[[:space:]]\" $variant_file>$tempdir/varfile";
        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
        );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }
        $variant_file = "$tempdir/varfile"
    }


    my %indelVariantHash;
    my %snvVariantHash;
    my %tooLongIndels;
    #store output variants so that we can check for missing ones at the end;
    my %foundHash;

    #read in all the variants and hash both the ref and var allele by position
    #also dump the snvs and indels in seperate files for readcounting
    my $inFh = IO::File->new( $variant_file ) || die "can't open file\n";
    open(SNVFILE,">$tempdir/snvpos");
    open(INDELFILE,">$tempdir/indelpos");
    while( my $sline = $inFh->getline )
    {
        chomp($sline);

        #skip header lines
        next if($sline =~ /^(#|Hugo_Symbol|Chr|chromosome)/i);

        my @fields = split("\t",$sline);

        $fields[3] =~ s/0/\-/;
        $fields[4] =~ s/0/\-/;
        $fields[3] =~ s/\*/\-/;
        $fields[4] =~ s/\*/\-/;
        $fields[3] = uc($fields[3]);
        $fields[4] = uc($fields[4]);


        my $key = join("\t",(@fields[0..1]));

        #is it an indel?
        if (($fields[3] =~ /\-/) || ($fields[4] =~ /\-/) ||
            (length($fields[3]) > 1) || (length($fields[4]) > 1)){

            #is it longer than the max length?
            if((length($fields[3]) > $indel_size_limit) || (length($fields[4]) > $indel_size_limit)){
                $tooLongIndels{join("\t",($fields[0],$fields[1],$fields[3],$fields[4]))} = 0;
            } else {
                #could have more than one indel per position
                if(defined($indelVariantHash{$key})){
                    $indelVariantHash{$key} = $indelVariantHash{$key} . "," . join("\t",($fields[3],$fields[4]));
                } else {
                    $indelVariantHash{$key} = join("\t",($fields[3],$fields[4]));
                }
                $foundHash{join("\t",($fields[0],$fields[1],$fields[3],$fields[4]))} = 0;
                print INDELFILE join("\t",($fields[0],$fields[1]-1,$fields[2],$fields[3],$fields[4])) . "\n";
            }

        } else { #snv
            #could have more than one snv per position
            if(defined($snvVariantHash{$key})){
                $snvVariantHash{$key} = $snvVariantHash{$key} . "," . join("\t",($fields[3],$fields[4]));
            } else {
                $snvVariantHash{$key} = join("\t",($fields[3],$fields[4]));
                print SNVFILE join("\t",($fields[0],$fields[1],$fields[2],$fields[3],$fields[4])) . "\n";
            }
            $foundHash{join("\t",($fields[0],$fields[1],$fields[3],$fields[4]))} = 0;
        }
    }
    close(INDELFILE);
    close(SNVFILE);
    close($inFh);



    #open the output file
    my $OUTFILE;
    open($OUTFILE,">$output_file") || die "can't open $output_file for writing\n";

    my @libraries;
    if($self->per_library) {
        my $bam_file = $self->bam_file;
        my @header = `samtools view -H $bam_file`;

        my %libraries = map { ($_) = $_ =~ /LB:(.+?)\s/; $_ => 1; } grep { /^\@RG/ } @header;
        @libraries = sort keys %libraries;

        print $OUTFILE "#chr\tpos\tref\tvar\t";
        for my $lib (@libraries) {
            print $OUTFILE join("\t", map { join("_", $lib, $_) } qw( ref_count var_count VAF )), "\t";
        }
        print $OUTFILE "\n";
    }



    #------------------------------------------
    #now run the readcounting on snvs
    if( -s "$tempdir/snvpos"){
        my $return = Genome::Model::Tools::Sam::Readcount->execute(
            use_version => 0.5,
            bam_file => $bam_file,
            minimum_mapping_quality => $min_mapping_quality,
            minimum_base_quality => $min_base_quality,
            output_file => "$tempdir/readcounts",
            reference_fasta => $fasta,
            region_list => "$tempdir/snvpos",
            per_library => $self->per_library,
        );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }

        #sort and dedup the bam-readcount output
        my $cmd_obj = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [ "$tempdir/readcounts" ],
            output_file => "$tempdir/readcounts.sorted",
            );
        $cmd_obj->execute;
        system( "uniq $tempdir/readcounts.sorted >$tempdir/readcounts.sorted.uniq" );

        #parse the results
        my $reader = new Genome::File::BamReadcount::Reader("$tempdir/readcounts.sorted.uniq");
        while(my $entry = $reader->next) {

            my $ref_count = 0;
            my $var_count = 0;
            my $var_freq = 0;
            my $knownRef;
            my $knownVar;

            my ($chr, $pos) = ($entry->chromosome, $entry->position);
            my $key = join("\t", $chr, $pos);
            if(!(defined($snvVariantHash{$key}))){
                print STDERR "WARNING: position $chr : $pos not found in input\n";
                next;
            }

            my @snvs = split(",",$snvVariantHash{$key});

            foreach my $pair (@snvs){
                my @as = split("\t",$pair);
                $knownRef = $as[0];
                $knownVar = $as[1];
                unless($self->per_library) {
                    my $var_freq = 0;

                    my ($ref_count, $var_count) = snvCounts($self,$entry, $knownRef, $knownVar);

                    if ($entry->depth ne '0') {
                        $var_freq = $var_count/$entry->depth * 100;
                    }

                    $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;

                    if($count_non_reference_reads){
                        $knownVar = "NonRef";
                    }

                    filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, $var_count, $var_freq,
                        $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                }
                else {
                    my %counts;
                    my $total_ref = 0;
                    my $total_var = 0;
                    my $total_var_freq = 0;

                    for my $lib ($entry->libraries) {
                        my $var_freq = 0;

                        my ($ref_count, $var_count) = $self->snvCounts($lib, $knownRef, $knownVar);

                        if ($lib->depth ne '0') {
                            $var_freq = $var_count/$lib->depth * 100;
                        }
                        $counts{$lib->name} = [$ref_count, $var_count, $var_freq];
                        $total_ref += $ref_count;
                        $total_var += $var_count;
                    }


                    $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;

                    if($count_non_reference_reads){
                        $knownVar = "NonRef";
                    }

                    if ($entry->depth ne '0') {
                        $total_var_freq = $total_var/$entry->depth * 100;
                    }
                    unless(shouldFilter($total_ref, $total_var, $total_var_freq, $min_depth, $max_depth, $min_vaf, $max_vaf)) {
                        my @ordered_counts;
                        for my $count_ref (@counts{@libraries}) {
                            if(defined $count_ref) {
                                push @ordered_counts, @$count_ref;
                            }
                            else {
                                push @ordered_counts, (0,0,0);
                            }
                        }
                        printLibs($OUTFILE, $chr, $pos, $knownRef, $knownVar, @ordered_counts);
                    }
                }
            }
        }
    }


    #--------------------------------------------
    #now indels, which gets tricky
    #the way pileup places the coordinates gets weird, so output the appropriate bases to look at:

    #if there are no indels, skip
    if( -s "$tempdir/indelpos"){
        my $return = Genome::Model::Tools::Sam::Readcount->execute(
            use_version => 0.5,
            bam_file => $bam_file,
            minimum_mapping_quality => $min_mapping_quality,
            minimum_base_quality => $min_base_quality,
            output_file => "$tempdir/readcounts_indel",
            reference_fasta => $fasta,
            region_list => "$tempdir/indelpos",
            insertion_centric => 1,
            per_library => $self->per_library,
        );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }

        #sort and dedup the bam-readcount output
        my $cmd_obj = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [ "$tempdir/readcounts_indel" ],
            output_file => "$tempdir/readcounts_indel.sorted",
            );
        $cmd_obj->execute;
        system( "uniq $tempdir/readcounts_indel.sorted >$tempdir/readcounts_indel.sorted.uniq" );

        #parse the results
        my $reader = new Genome::File::BamReadcount::Reader("$tempdir/readcounts_indel.sorted.uniq");
        while( my $entry = $reader->next) {

            my $key = join("\t", $entry->chromosome, $entry->position);
            unless((defined($indelVariantHash{$key}))){
                next;
            }

            #can be more than one allele at each position
            my @alleles = split(",",$indelVariantHash{$key});
            foreach my $allele (@alleles){
                my ($knownRef, $knownVar) = split("\t",$allele);
                my $testvarallele;
                if($knownRef =~ /0|\-|\*/) { #INS
                    $testvarallele = "+$knownVar"
                }
                elsif ($knownVar =~ /0|\-|\*/){ #DEL
                    $testvarallele = "-$knownRef"
                }
                else {
                    print "WARNING - $knownRef/$knownVar isn't an indel, how did it get in the indel hash?\n";
                }

                unless($self->per_library) {
                    my $ref_count = 0;
                    my $var_count = 0;
                    ($ref_count, $var_count) = $self->indelCounts($entry, $testvarallele);
                    if ($entry->depth ne '0') {
                        filterAndPrint($entry->chromosome, $entry->position, $knownRef, $knownVar, $ref_count, $var_count, ($var_count/$entry->depth)*100,
                                       $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                    }
                    else {
                        filterAndPrint($entry->chromosome, $entry->position, $knownRef, $knownVar, $ref_count, $var_count, 0,
                                       $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                    }
                    $foundHash{join("\t",$entry->chromosome,$entry->position,$knownRef,$knownVar)} = 1;
                }
                else {
                    my %counts;
                    my $total_ref = 0;
                    my $total_var = 0;
                    my $total_var_freq = 0;
                    for my $lib ($entry->libraries) {
                        my $var_freq = 0;
                        my ($ref_count, $var_count) = $self->indelCounts($lib, $testvarallele);
                        if ($lib->depth ne '0') {
                            $var_freq = $var_count/$lib->depth * 100;
                        }
                        $counts{$lib->name} = [$ref_count, $var_count, $var_freq];
                        $total_ref += $ref_count;
                        $total_var += $var_count;
                    }
                    $foundHash{join("\t",$entry->chromosome,$entry->position,$knownRef,$knownVar)} = 1;
                    if ($entry->depth ne '0') {
                        $total_var_freq = $total_var/$entry->depth * 100;
                    }
                    unless(shouldFilter($total_ref, $total_var, $total_var_freq, $min_depth, $max_depth, $min_vaf, $max_vaf)) {
                        #this is duplicated from above and should be refactored
                        my @ordered_counts;
                        for my $count_ref (@counts{@libraries}) {
                            if(defined $count_ref) {
                                push @ordered_counts, @$count_ref;
                            }
                            else {
                                push @ordered_counts, (0,0,0);
                            }
                        }
                        printLibs($OUTFILE, $entry->chromosome, $entry->position, $knownRef, $knownVar, @ordered_counts);
                    }
                }
            }
        }
    }

    #Check all the variants and output those with no output (had 0 reads)
    foreach my $k (keys(%foundHash)){
        unless($foundHash{$k}){
            #site not called, gets a zero count
            my ($chr, $pos, $knownRef, $knownVar) = split("\t",$k);
            if($self->per_library) {
                my $num_libs = scalar(@libraries) || 1;
                printLibs($OUTFILE, $chr, $pos, $knownRef, $knownVar, (0) x ($num_libs * 3));
            }
            else {
                filterAndPrint($chr, $pos, $knownRef, $knownVar, 0, 0, 0,
                    $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
            }
        }
    }

    foreach my $k (keys(%tooLongIndels)){
        #site too long, gets an NA value
        my ($chr, $pos, $knownRef, $knownVar) = split("\t",$k);
        if($self->per_library) {
            my $num_libs = scalar(@libraries) || 1;
            printLibs($OUTFILE, $chr, $pos, $knownVar, $knownVar, ("NA") x ($num_libs * 3));
        }
        else {
            filterAndPrint($chr, $pos, $knownRef, $knownVar, "NA", "NA", "NA",
                $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
        }
    }
    close($OUTFILE);

    return(1);
}

sub indelCounts {
    my ($self, $lib, $testvarallele) = @_;
    my $ref_count = 0;
    my $var_count = 0;
    for my $allele ($lib->alleles) {
        if($allele ne $testvarallele) {
            $ref_count += $lib->metrics_for($allele)->count;
        }
        else {
            $var_count += $lib->metrics_for($allele)->count;
        }
    }
    return ($ref_count, $var_count);
}

sub snvCounts {
    my ($self, $lib, $knownRef, $knownVar) = @_;
    my $ref_count = 0;
    my $var_count = 0;

    for my $allele ($lib->alleles) {
        # assume that the ref call is ACTG, not iub
        # (assumption looks valid in my files)
        if ($allele eq $knownRef){
            $ref_count += $lib->metrics_for($allele)->count;
            next;
        }

        # if we're counting all non-reference reads, not just the specified allele
        if($self->count_non_reference_reads){
            unless($allele eq $knownRef){
                $var_count += $lib->metrics_for($allele)->count;
            }
            next;
        }

        # if this base is included in the IUB code for
        # for the variant, (but doesn't match the ref)
        if (matchIub($allele,$knownRef,$knownVar)){
            $var_count += $lib->metrics_for($allele)->count;
        }
    }
    return ($ref_count, $var_count);
}
