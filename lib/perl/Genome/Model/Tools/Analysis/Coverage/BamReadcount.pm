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

    use_varscan => {
        is => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => 'use samtools mpilup and varscan readcounts for snv readcounts'
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
    my $use_varscan = $self->use_varscan;

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
            return 1;
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
                if(length($fields[3]) > 1 && !$use_varscan) {
                    $fields[1] -= 1;
                }
                print INDELFILE join("\t",($fields[0],$fields[1],$fields[2],$fields[3],$fields[4])) . "\n";
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

        # use samtools, then varscan for gathering readcounts #
        if ($use_varscan) {

            open(SNVMOD,">$tempdir/snvpos.bed");
            $inFh = IO::File->new( "$tempdir/snvpos" ) || die "can't open file\n";

            #convert the snv file to bed
            while( my $line = $inFh->getline )
            {
                chomp($line);
                my @F = split("\t",$line);
                $F[1]--;
                print SNVMOD join("\t",@F) . "\n";
            }
            close(SNVMOD);

            #run samtools view, then mpileup to get the readcounts for each snv:
            my $cmd = "samtools view -ub -L $tempdir/snvpos.bed $bam_file | samtools mpileup -f $fasta -q $min_mapping_quality - > $tempdir/snv_mpileup";
            $self->status_message("Running command: $cmd");

            my $return = Genome::Sys->shellcmd(
                cmd => "$cmd",
            );
            unless($return) {
                $self->error_message("Failed to execute: Returned $return");
                die $self->error_message;
            }

            #run varscan to parse the samtools file
            $cmd = "java -jar /gsc/scripts/lib/java/VarScan/VarScan.v2.2.9.jar readcounts $tempdir/snv_mpileup --min-coverage 1 --min-base-qual $min_base_quality --output-file $tempdir/snvs.varscan";
            $self->status_message("Running command: $cmd");

            $return = Genome::Sys->shellcmd(
                cmd => "$cmd",
            );
            unless($return) {
                $self->error_message("Failed to execute: Returned $return");
                die $self->error_message;
            }

            #read, clean up and print snvs
            $inFh = IO::File->new( "$tempdir/snvs.varscan" ) || die "can't open varscan snv file\n";
            while( my $line = $inFh->getline )
            {
                chomp($line);
                next if $line =~ /^chrom/;
                my ($chr, $pos, $ref, $depth, $q0_depth, @counts) = split("\t",$line); #counts format: base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads

                my $ref_count = 0;
                my $var_count = 0;
                my $var_freq = 0;
                my $knownRef;
                my $knownVar;

                if(!(defined($snvVariantHash{join("\t",($chr, $pos))}))){
                    #print STDERR "WARNING: position $chr : $pos not found in input\n";
                    next;
                }

                my @snvs = split(",",$snvVariantHash{join("\t",($chr, $pos))});

                foreach my $pair (@snvs){
                    my @as = split("\t",$pair);
                    $knownRef = $as[0];
                    $knownVar = $as[1];
                    my $ref_count = 0;
                    my $var_count = 0;
                    my $var_freq = 0;

                    #go through each base at that position, grab the correct one
                    foreach my $count_stats (@counts) {
                        my ($allele, $count, $strands, $bq, $mq, $plus_reads, $minus_reads) = split /:/, $count_stats;

                        # assume that the ref call is ACTG, not iub
                        # (assumption looks valid in my files)
                        if ($allele eq $knownRef){
                            $ref_count += $count;
                            next;
                        }

                        # if we're counting all non-reference reads, not just the specified allele
                        if($count_non_reference_reads){
                            unless($allele eq $knownRef){
                                $var_count += $count;
                            }
                            next;
                        }

                        # if this base is included in the IUB code for
                        # for the variant, (but doesn't match the ref)
                        if (matchIub($allele,$knownRef,$knownVar)){
                            $var_count += $count;
                        }

                    }
                    if ($depth ne '0') {
                        $var_freq = $var_count/$depth * 100;
                    }

                    $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;

                    if($count_non_reference_reads){
                        $knownVar = "NonRef";
                    }

                    filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, $var_count, $var_freq,
                        $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);

                }
            }
        }

        # else use bam-readcount #
        else {
            my $return = Genome::Model::Tools::Sam::Readcount->execute(
                use_version => 0.5,
                bam_file => $bam_file,
                minimum_mapping_quality => $min_mapping_quality,
                minimum_base_quality => $min_base_quality,
                output_file => "$tempdir/readcounts",
                reference_fasta => $fasta,
                region_list => "$tempdir/snvpos",
            );
            unless($return) {
                $self->error_message("Failed to execute: Returned $return");
                die $self->error_message;
            }

            #parse the results
            my $reader = new Genome::File::BamReadcount::Reader("$tempdir/readcounts");
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
                    my $ref_count = 0;
                    my $var_count = 0;
                    my $var_freq = 0;

                    #go through each base at that position, grab the correct one
                    for my $allele ($entry->alleles) {
                        # assume that the ref call is ACTG, not iub
                        # (assumption looks valid in my files)
                        if ($allele eq $knownRef){
                            $ref_count += $entry->metrics_for($allele)->count;
                            next;
                        }

                        # if we're counting all non-reference reads, not just the specified allele
                        if($count_non_reference_reads){
                            unless($allele eq $knownRef){
                                $var_count += $entry->metrics_for($allele)->count;
                            }
                            next;
                        }

                        # if this base is included in the IUB code for
                        # for the variant, (but doesn't match the ref)
                        if (matchIub($allele,$knownRef,$knownVar)){
                            $var_count += $entry->metrics_for($allele)->count;
                        }

                    }
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
            }
        }
    }


    #--------------------------------------------
    #now indels, which gets tricky
    #the way pileup places the coordinates gets weird, so output the appropriate bases to look at:

    #if there are no indels, skip
    if( -s "$tempdir/indelpos"){
        if($self->use_varscan) {

            #grab pileups for each indel
            $inFh = IO::File->new( "$tempdir/indelpos" ) || die "can't open file\n";
            while( my $line = $inFh->getline )
            {
                #convert coordinates to grab the appropriate bases:
                chomp($line);
                my @F = split("\t",$line);
                if($F[3] =~ /0|\-|\*/ ){ #INS
                    $F[1]--;
                    $F[2]--;
                } elsif ($F[4] =~ /0|\-|\*/){ #DEL - get two bases, since the del is placed on the previous base
                    $F[1] = $F[1] - 2;
                    $F[2] = $F[1] + 2;
                } else {
                    print STDERR "WARNING: bad indel format: $line \n";
                }

                #run mpileup to get the readcounts:
                my $cmd = "samtools mpileup -f $fasta -q 1 -r $F[0]:$F[1]-$F[2]  $bam_file >>$tempdir/pileup";
                #$self->status_message("Running command: $cmd");

                my $return = Genome::Sys->shellcmd(
                    cmd => "$cmd",
                );
                unless($return) {
                    $self->error_message("mpileup failure. Tried to run:\n$cmd");
                    $self->error_message("Returned:\n$return");
                    die $self->error_message;
                }
            }

            #run varscan to parse the samtools file
            my $cmd = "java -jar /gsc/scripts/lib/java/VarScan/VarScan.v2.2.9.jar readcounts $tempdir/pileup --min-coverage 1 --min-base-qual $min_base_quality --output-file $tempdir/indels.varscan";
            $self->status_message("Running command: $cmd");

            my $return = Genome::Sys->shellcmd(
                cmd => "$cmd",
            );
            unless($return) {
                $self->error_message("Failed to execute: Returned $return");
                die $self->error_message;
            }

            my %readdepth;
            my %reads;

            #store indel counts
            $inFh = IO::File->new( "$tempdir/indels.varscan" ) || die "can't open varscan file\n";
            while( my $line = $inFh->getline )
            {
                chomp($line);
                next if $line =~ /^chrom/;
                my @F =  split("\t",$line);

                my $chr = $F[0];
                my $pos = $F[1];
                my $key = join("\t",($chr,$pos));

                $readdepth{$key} = $F[3];
                for my $i (5..$#F){
                    if($F[$i] =~ /\:/){ #skip blanks
                        my @counts = split(":",$F[$i]);
                        my $base = $counts[0];
                        $reads{$key}{$base} = $counts[1];
                    }
                }
            }

            #go through the indels we're looking for and grab their depths
            foreach my $key (keys(%indelVariantHash)){
                my ($chr, $pos) = split("\t",$key);
                my ($refbase,$varbase) = split("\t",$indelVariantHash{$key});

                if($refbase =~ /0|\-|\*/ ){ #INS
                    #insertion is easier, just take the insertion count divided by depth for vaf
                    my $depth = 0;
                    my $refcount = 0;
                    my $varcount = 0;

                    if(defined($readdepth{$key})){
                        $depth = $readdepth{$key};

                        foreach my $base (keys(%{$reads{$key}})){
                            if($base =~ /INS-\d+-$varbase/){
                                $varcount = $reads{$key}{$base};
                            }
                        }

                        #all reads will include ref for an insertion in mpileup-land, so
                        #we subtract the var from the ref
                        $refcount = $refcount - $varcount;
                        filterAndPrint($chr, $pos, $refbase, $varbase, $depth-$varcount, $varcount, ($varcount/$depth)*100,
                            $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                    } else {
                        #if it wasn't in the pileup, it wasn't covered, so vals will remain zero.
                        filterAndPrint($chr, $pos, $refbase, $varbase, $refcount, $varcount, 0,
                            $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                    }


                } elsif ($varbase =~ /0|\-|\*/){ #DEL
                    #deletions are tricky. The deletion gets placed on the previous base, but the ref counts 
                    #and depth are on the correct base.
                    my $testrefbase = $refbase;
                    if (length($refbase) > 1){
                        $testrefbase = substr($refbase,0,1);
                    }

                    my $depth = 0;
                    my $refcount = 0;
                    my $varcount = 0;
                    if(defined($readdepth{$key})){
                        my $depth = $readdepth{$key};

                        #first check the correct base
                        foreach my $base (keys(%{$reads{$key}})){
                            if($base eq $testrefbase){
                                $refcount = $reads{$key}{$base};
                            }
                        }
                        #now check the preceding base
                        my $pkey = join("\t",($chr,$pos-1));
                        foreach my $base (keys(%{$reads{$pkey}})){
                            if($base =~ /DEL-\d+-$testrefbase/){
                                $varcount = $reads{$pkey}{$base};
                            }
                        }

                        filterAndPrint($chr, $pos, $refbase, $varbase, $refcount, $varcount, ($varcount/$depth)*100,
                            $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);            
                    } else {
                        #if it wasn't in the pileup, it wasn't covered, or didn't exist, 
                        #so ref base goes to depth
                        filterAndPrint($chr, $pos, $refbase, $varbase, $refcount, $varcount, 0,
                            $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                    }                
                } else {
                    print "WARNING - $refbase/$varbase isn't an indel, how did it get in the indel hash?\n";
                }

                $foundHash{join("\t",$chr,$pos,$refbase,$varbase)} = 1;

            }
        }
        else {

            my $return = Genome::Model::Tools::Sam::Readcount->execute(
                use_version => 0.5,
                bam_file => $bam_file,
                minimum_mapping_quality => $min_mapping_quality,
                minimum_base_quality => $min_base_quality,
                output_file => "$tempdir/readcounts_indel",
                reference_fasta => $fasta,
                region_list => "$tempdir/indelpos",
                insertion_centric => 1,
            );
            unless($return) {
                $self->error_message("Failed to execute: Returned $return");
                die $self->error_message;
            }

            my $reader = new Genome::File::BamReadcount::Reader("$tempdir/readcounts_indel");
            while( my $entry = $reader->next) { 

                my $ref_count = 0;
                my $var_count = 0;

                my $key = join("\t", $entry->chromosome, $entry->position);
                unless((defined($indelVariantHash{$key}))){
                    next;
                }
                my ($knownRef, $knownVar) = split("\t",$indelVariantHash{$key});
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

                for my $allele ($entry->alleles) {
                    if($allele ne $testvarallele) {
                        $ref_count += $entry->metrics_for($allele)->count;
                    }
                    else {
                        $var_count += $entry->metrics_for($allele)->count;
                    }
                }
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
        }
    }
    
    #Check all the variants and output those with no output (had 0 reads)
    foreach my $k (keys(%foundHash)){
        unless($foundHash{$k}){
            #site not called, gets a zero count
            my ($chr, $pos, $knownRef, $knownVar) = split("\t",$k);
            filterAndPrint($chr, $pos, $knownRef, $knownVar, 0, 0, 0,
                           $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
        }
    }

    foreach my $k (keys(%tooLongIndels)){
        #site too long, gets an NA value
        my ($chr, $pos, $knownRef, $knownVar) = split("\t",$k);
        filterAndPrint($chr, $pos, $knownRef, $knownVar, "NA", "NA", "NA",
                       $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
    }
    close($OUTFILE);
    
    return(1);
}

