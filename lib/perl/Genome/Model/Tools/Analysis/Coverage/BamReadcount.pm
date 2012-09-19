package Genome::Model::Tools::Analysis::Coverage::BamReadcount;
use strict;
use Genome;
use IO::File;
use warnings;


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
            default => 2,
        }

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
    my $min_quality_score = $self->min_quality_score;

    my $min_vaf = $self->min_vaf;
    my $max_vaf = $self->max_vaf;
    my $min_depth = $self->min_depth;
    my $max_depth = $self->max_depth;

    my $chrom = $self->chrom;

    my $fasta;
    if ($genome_build eq "36") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->data_directory . "/all_sequences.fa";
    }
    elsif (($genome_build eq "37") || ($genome_build eq "37lite")) {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $fasta = $reference_build_fasta_object->data_directory . "/all_sequences.fa";
    }
    elsif ($genome_build eq "mus37") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-mouse-build37");
        $fasta = $reference_build_fasta_object->data_directory . "/all_sequences.fa";        
    } elsif ($genome_build eq "mus37wOSK") {
        $fasta = "/gscmnt/sata135/info/medseq/dlarson/iPS_analysis/lentiviral_reference/mousebuild37_plus_lentivirus.fa";
    } elsif ($genome_build eq "mm9") {
        $fasta = "/gscmnt/ams1102/info/model_data/2869962191/build107494762/all_sequences.fa";
    } elsif (-e $genome_build ) {
        $fasta = $genome_build;
    } else {
        die ("invalid genome build or fasta path: $genome_build\n");
    }



    #--------------------------------------------------
    #convert iub bases to lists
    sub convertIub{
        my ($base) = @_;
        
        #deal with cases like "A/T" or "C/W"
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

        # use a lookup table to return the correct base
        # there's a more efficient way than defining this, 
        # every time, but meh.
        my %iub_codes;
        $iub_codes{"A"}="A";
        $iub_codes{"C"}="C";
        $iub_codes{"G"}="G";
        $iub_codes{"T"}="T";
        $iub_codes{"U"}="T";
        $iub_codes{"M"}="A,C";
        $iub_codes{"R"}="A,G";
        $iub_codes{"W"}="A,T";
        $iub_codes{"S"}="C,G";
        $iub_codes{"Y"}="C,T";
        $iub_codes{"K"}="G,T";
        $iub_codes{"V"}="A,C,G";
        $iub_codes{"H"}="A,C,T";
        $iub_codes{"D"}="A,G,T";
        $iub_codes{"B"}="C,G,T";
        $iub_codes{"N"}="A,C,G,T";

        return $iub_codes{$base}
    }


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
    
    
    sub filterAndPrint{
        my ($chr, $pos, $knownRef, $knownVar, $ref_count, $var_count, $var_freq, 
         $min_depth,$max_depth,$min_vaf,$max_vaf,$OUTFILE) = @_;
        #filters on the output
        my $do_print=1;

        #handle the special case where this value is NA, which means don't filter, but
        #pass everything through. This lets AddReadcounts.pm work correctly.
        if(defined($min_depth)){
            unless($min_depth eq "NA"){
                if($var_freq eq "NA"){
                    $do_print = 0;
                } else {
                    $do_print = 0 if(($ref_count + $var_count) < $min_depth);
                }
            }
        }                
        if(defined($max_depth)){
            unless($max_depth eq "NA"){
                if($var_freq eq "NA"){
                    $do_print = 0;
                } else {
                    $do_print = 0 if(($ref_count + $var_count) > $max_depth);
                }
            }
        }                
 
        if(defined($min_vaf)){            
            unless($min_vaf eq "NA"){
                if($var_freq eq "NA"){
                    $do_print = 0;
                } elsif( $var_freq < $min_vaf) {
                    $do_print = 0;
                }
            }
        }
        if(defined($max_vaf)){            
            unless($max_vaf eq "NA"){
                if($var_freq eq "NA"){
                    $do_print = 0;
                } elsif( $var_freq > $max_vaf) {
                    $do_print = 0;
                }
            }
        }

        
        if($do_print){
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
            if((length($fields[3]) > $self->indel_size_limit) || (length($fields[4]) > $self->indel_size_limit)){
                $tooLongIndels{join("\t",($fields[0],$fields[1],$fields[3],$fields[4]))} = 0;
            } else {
                #could have more than one indel per position
                if(defined($indelVariantHash{$key})){
                    $indelVariantHash{$key} = $indelVariantHash{$key} . "," . join("\t",($fields[3],$fields[4]));
                } else {
                    $indelVariantHash{$key} = join("\t",($fields[3],$fields[4]));
                }
                $foundHash{join("\t",($fields[0],$fields[1],$fields[3],$fields[4]))} = 0;
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
    
  
    #------------------------------------------
    #now run the readcounting on snvs
    if( -s "$tempdir/snvpos"){
        my $cmd = "bam-readcount -q $min_quality_score -f $fasta -l $tempdir/snvpos $bam_file >$tempdir/readcounts";
        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }

        #parse the results
        my $inFh2 = IO::File->new( "$tempdir/readcounts" ) || die "can't open file\n";
        while( my $line = $inFh2->getline )
        {
            chomp($line);
            my ($chr, $pos, $ref, $depth, @counts) = split("\t",$line);
            
            my $ref_count = 0;
            my $var_count = 0;
            my $var_freq = 0;
            my $knownRef;
            my $knownVar;
            
            if(!(defined($snvVariantHash{join("\t",($chr, $pos))}))){
                print STDERR "WARNING: position $chr : $pos not found in input\n";
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
                    my ($allele, $count, $mq, $bq) = split /:/, $count_stats;
                    
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
                }                        
                filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, $var_count, $var_freq,
                               $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE);
                
                $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;
            }
            
        }
    }    
    
    #--------------------------------------------
    #now indels, which gets tricky
    #the way pileup places the coordinates gets weird, so output the appropriate bases to look at:

    #if there are no indels, this could happen
    if( -s "$tempdir/indelpos"){
        open(INDELMOD,">$tempdir/indelpos.mod");
        $inFh = IO::File->new( "$tempdir/indelpos" ) || die "can't open file\n";

        #convert the indel file to bed
        while( my $line = $inFh->getline )
        {
            chomp($line);
            my @F = split("\t",$line);
            if($F[3] =~ /0|\-|\*/ ){ #INS
                $F[1]--;
                $F[2]--;
            } elsif ($F[4] =~ /0|\-|\*/){ #DEL
                $F[1] = $F[1] - 2;
                $F[2] = $F[1] + 1;
            } else {
                print STDERR "WARNING: bad indel format: $line \n";
            }

            print INDELMOD join("\t",@F) . "\n";
        }
        close(INDELMOD);

        #run mpileup to get the readcounts for each indel:
        my $cmd = "samtools mpileup -f $fasta -q 1 -l $tempdir/indelpos.mod $bam_file >$tempdir/pileup";
        $self->status_message("Running command: $cmd");

        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }

        #run varscan to parse the samtools file
        $cmd = "java -jar /gsc/scripts/lib/java/VarScan/VarScan.v2.2.9.jar readcounts $tempdir/pileup --min-coverage 1 --min-base-qual 0 --output-file $tempdir/indels.varscan";
        $self->status_message("Running command: $cmd");
        
        $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }       
        
        #read, clean up and print indels
        $inFh = IO::File->new( "$tempdir/indels.varscan" ) || die "can't open varscan file\n";
        while( my $line = $inFh->getline )
        {
            chomp($line);
            next if $line =~ /^chrom/;
            my @F =  split("\t",$line);
            
            my $ref_count = 0;
            my $var_count = 0;
            my $var_freq = 0;
            my $knownRef;
            my $knownVar;
            my $chr = $F[0];
            my $pos = $F[1];
            my @indels;

            if(defined($indelVariantHash{join("\t",($chr, $pos))})){
                @indels = split(",",$indelVariantHash{join("\t",($chr, $pos))}); 
            } elsif (defined($indelVariantHash{join("\t",($chr, $pos+1))})){
                $pos++;
                @indels = split(",",$indelVariantHash{join("\t",($chr, $pos))});
            } else {
                print STDERR "WARNING: position $chr : $pos not found in input\n";
                next;
            }

            #can be more than one indel per position

            foreach my $pair (@indels){
                my @as = split("\t",$pair);
                $knownRef = $as[0];
                $knownVar = $as[1];                

                my @refdetails = split(":",$F[5]);
                $ref_count = $refdetails[1];
                
                #can be more than one variant output in trailing columns of varscan format
                my $found = 0;
                my $i;
                for($i=6;$i<@F;$i++){
                    my @varResults = split(":",$F[$i]);
                    my $varReads = $varResults[1];
                    
                    if ($varResults[0] =~ /INS/){
                        my @varInfo = split("-",$varResults[0]);
                        my $varbase = $varInfo[2];
                        if($varbase eq $knownVar){
                            unless ($varResults[1] == 0){
                                $var_freq = ($varResults[1]/$F[3])*100;
                            }
                            filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, $varResults[1], $var_freq,
                                           $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE); 
                            
                            $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;
                            $found = 1;
                        }                    
                    } elsif ($varResults[0] =~ /DEL/){
                        my @varInfo = split("-",$varResults[0]);
                        my $varbase = $varInfo[2];
                        if($varbase eq $knownRef){
                            unless ($varResults[1] == 0){
                                #this isn't entirely correct. To get the proper depth, we should also
                                #query the following base, get the depth and use it (since 
                                #samtools reports deletions on the previous base). For now, we're
                                #using the depth at the prev base as a proxy. Should be accurate to
                                #within a read or two.
                                $var_freq = ($varResults[1]/($F[3]+$varResults[1]))*100;
                            }
                            filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, $varResults[1], $var_freq,
                                           $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE); 
                            $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;
                            $found = 1;
                        }
                    }
                }
                #if the site is in the output, but the indel isn't there, we output zero for the var count
                unless($found){
                    filterAndPrint($chr, $pos, $knownRef, $knownVar, $ref_count, 0, 0,
                                   $min_depth, $max_depth, $min_vaf, $max_vaf, $OUTFILE); 
                    $foundHash{join("\t",$chr,$pos,$knownRef,$knownVar)} = 1;
                }
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
}
