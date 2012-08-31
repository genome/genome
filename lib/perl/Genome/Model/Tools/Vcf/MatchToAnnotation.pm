package Genome::Model::Tools::Vcf::MatchToAnnotation;

use strict;
use warnings;
use FileHandle;
use Genome;

class Genome::Model::Tools::Vcf::MatchToAnnotation {
    is => 'Command',

    has => [
        vcf_file   => {
            is => 'Text',
            doc => "VCF file in uncompressed format. Assumes there are two samples included (tumor/normal)",
            is_input => 1},

        annotation_file => {
            is => 'Text',
            doc => "WU annotation file containing variants that passed validation (at least 5 col, 1 based)",
            is_input => 1},

        output_file   => {
            is => 'Text',
            doc => "VCF file with filters updated",
            is_input => 1,
            is_output => 1},

        somatic_variation_model_id => {
            is => 'Text',
            doc => "if the annotation file has variants that aren't in the VCF, provide a somVar model id so that we can backfill those variant locations",
            is_input => 1,
            is_output => 1,
            optional => 1},

        ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Update a VCF file to match validated variants in an annotation file"
}

sub help_synopsis {
    return <<EOS
        This command takes a VCF output by the pipeline, along with an annotation file containing validated somatic variants. It updates the filter field to make sure that everything listed as passing in the anno file is marked PASS and somatic, and anything not listed in the validated file is marked as non-passing.
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS
        This command takes a VCF output by the pipeline, along with an annotation file containing validated somatic variants. It updates the filter field to make sure that everything listed as passing in the anno file is marked PASS and somatic, and anything not listed in the validated file is marked as non-passing.
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $vcf_file = $self->vcf_file;
    my $anno_file = $self->annotation_file;
    my $output_file = $self->output_file;
    my $somvar_model_id = $self->somatic_variation_model_id;


    ## Load the annotation file ##
    my %annohash;
    my $inFh = IO::File->new( $anno_file ) || die "can't open anno file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        my $key = join("\t",(@F[0..1],@F[3..4]));
        $annohash{$key} = 0;
    }
    close($inFh);


    open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

    if(defined($somvar_model_id)){
        open(TEMPFILETUMOR, ">$output_file.tumor.tofill") or die "Can't open outfile: $!\n";
        open(TEMPFILENORMAL, ">$output_file.normal.tofill") or die "Can't open outfile: $!\n";
    }


    #schlep through the VCF file
    my $linecount = 0;
    my $firstline;
    my $ft_pos = -1;
    $inFh = IO::File->new( $vcf_file ) || die "can't open anno file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);

        my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field, $field2) = split(/\t/, $line);

        #print the header
        if ($line =~ /^#/) {
            print OUTFILE $line . "\n";
            if(defined($somvar_model_id)){
                if($line =~ /^#CHROM/){
                    print TEMPFILETUMOR join("\t",($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field)) . "\n";
                    print TEMPFILENORMAL join("\t",($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field2)) . "\n"; 
                } else {
                    print TEMPFILETUMOR $line . "\n";
                    print TEMPFILENORMAL $line . "\n";
                }
            } 
            next;
        }

        #figure out where the FT field is
        if($linecount == 0){
            my @formats = split(":",$format);
            for(my $i=0; $i<@formats; $i++){
                if($formats[$i] eq "FT"){
                    $ft_pos = $i;
                    last;
                }
            }
            $firstline = $line;
        }

        my @fields = split(":",$field);
        # my @fields2 = split(":",$field2);

        #if ft field doesn't exist, add it to the end of the fields/format lines
        if($ft_pos == -1){
            $format = $format . ":FT";
            $ft_pos = @fields
        }
        $linecount++;



        my @alts = split(",",$var);
        foreach my $a (@alts){
            if(defined($annohash{join("\t",($chrom, $position, $ref, $a))})){
                #valid, make sure it's marked pass
                unless ($fields[$ft_pos] =~ /PASS/){
                    $fields[$ft_pos] = "PASS";
                }
                # unless($field2 eq "."){
                #     unless ($fields2[$ft_pos] =~ /PASS/){
                #         $fields2[$ft_pos] = "PASS";
                #     }
                # }
                #mark as found
                $annohash{join("\t",($chrom, $position, $ref, $a))} = 1;
            } else {
                #invalid, make sure it's marked fail
                if ($fields[$ft_pos] =~ /PASS/){
                    $fields[$ft_pos] = "DidNotValidate";
                }
                # unless($field2 eq "."){
                #     if ($fields2[$ft_pos] =~ /PASS/){
                #         $fields2[$ft_pos] = "DidNotValidate";
                #     }
                # }
            }
        }

        $field = join(":",@fields);
        # $field2 = join(":",@fields2);


        print OUTFILE join("\t",$chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field, $field2) . "\n";
    }
    close($inFh);


    # Add a third category - sites missing from the vcf. 
    # if a som-var model is defined, split these into single-sample VCFs 
    # so gmt vcf backfill will work. Otherwise, just fill with dots.
    my $added = 0;
    foreach my $k (keys(%annohash)){
        unless($annohash{$k}){
            #use the first line as a template
            my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field, $field2) = split(/\t/, $firstline);
            my ($kchrom, $kposition, $kref, $kvar) = split("\t",$k);

            #don't have quality info, so we just fill it with dots
            $chrom = $kchrom;
            $position = $kposition;
            $id = ".";
            $ref = $kref;
            $var = $kvar;
            $score = ".";
            $filter = ".";
            $info = ".";

            my @fields = split(":",$field);
            for(my $i=0; $i<@fields; $i++){
                if($i == $ft_pos){
                    $fields[$i] = "PASS";
                } else {
                    my @vals = split("/",$fields[$i]);
                    for(my $j=0; $j<@vals; $j++){
                        $vals[$j] = ".";
                    }
                    $fields[$i] = join("/",@vals);
                }
            }
            $field = join(":",@fields);
                        
            unless($field2 eq "."){
                my @fields2 = split(":",$field2);
                for(my $i=0; $i<@fields2; $i++){
                    if($i == $ft_pos){
                        $fields2[$i] = "PASS";
                    } else {
                        my @vals = split("/",$fields2[$i]);
                        for(my $j=0; $j<@vals; $j++){
                            $vals[$j] = ".";
                        }
                        $fields2[$i] = join("/",@vals);
                    }
                }
                $field2 = join(":",@fields2);
            }

            if(defined($somvar_model_id)){
                print TEMPFILETUMOR join("\t",($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field)) . "\n";
                print TEMPFILENORMAL join("\t",($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field2)) . "\n";
            } else {
                print OUTFILE join("\t",$chrom, $position, $id, $ref, $var, $score, $filter, $info, $format, $field, $field2) . "\n";
            }
            
            $added = 1;
        }
    }
    
    close(TEMPFILETUMOR);
    close(TEMPFILENORMAL);

    # backfill VCFs if we've added anything
    if($added && defined($somvar_model_id))
    {
        my $model = Genome::Model->get($somvar_model_id) or die "Could not find model ($somvar_model_id\n";

        my $tumor_model = $model->tumor_model;
        my $normal_model = $model->normal_model;

        my $tumor_build = $tumor_model->last_succeeded_build;
        my $normal_build = $normal_model->last_succeeded_build;

        
        my $tumor_samtools = glob($tumor_build->data_directory . "/variants/snv/samtools*/snvs.hq");
        my $normal_samtools = glob($normal_build->data_directory . "/variants/snv/samtools*/snvs.hq");

        
        print STDERR "tumor: $tumor_samtools\n";
        print STDERR "normal: $normal_samtools\n";
        
        #sort the output
        `vcf-sort $output_file.tumor.tofill>$output_file.tumor.tofill.sorted`;
        `vcf-sort $output_file.normal.tofill>$output_file.normal.tofill.sorted`;
        `mv -f $output_file.tumor.tofill.sorted $output_file.tumor.tofill`;
        `mv -f $output_file.normal.tofill.sorted $output_file.normal.tofill`;
        #get the list of variants
        `grep -v "^#" $output_file.tumor.tofill | cut -f 1,2,5 >$output_file.posfile`;

        # my $cmd = Genome::Model::Tools::Vcf::Backfill->create(
        #     bam_file => $bam_file,
        #     output_file =>  "$tempdir/rcfile",
        #     variant_file => $variant_file,
        #     genome_build => $genome_build, 
        #     chrom => $chrom,
        #     min_depth  => $min_depth,
        #     max_depth => $max_depth,
        #     min_vaf => $min_vaf,
        #     max_vaf => $max_vaf,
        #     indel_size_limit => $indel_size_limit,
        # );
        # unless ($cmd->execute) {
        #     die "Bam-readcount failed";
        # }
        # gmt vcf backfill
        
    }      


    #sort the output with joinx if we've added anything
    if($added){
        `vcf-sort $output_file >$output_file.sorted && mv -f $output_file.sorted $output_file`;
    }
    return 1;
}
