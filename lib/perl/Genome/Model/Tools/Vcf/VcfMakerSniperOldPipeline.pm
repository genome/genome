package Genome::Model::Tools::Vcf::VcfMakerSniperOldPipeline;

use strict;
use warnings;
use Genome;
use File::stat;
#use Time::localtime;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
#use POSIX qw(log10);
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);
use Data::Dumper;
use Genome::Model::Tools::Vcf::Helpers qw/genGT/;


class Genome::Model::Tools::Vcf::VcfMakerSniperOldPipeline {
    is => 'Command',
    has => [
    output_file => {
        is => 'Text',
        is_output => 1,
        doc => "List of mutations in Vcf format",
    },

    tumor_bam_file => {
        is => 'Text',
        doc => "Tumor sample bam file (for header)" ,
        is_optional => 0,
        is_input => 1},

    normal_bam_file => {
        is => 'Text',
        doc => "Normal sample bam file (for header)" ,
        is_optional => 0,
        is_input => 1},

    tumor_snp_file => {
        is => 'Text',
        doc => "all snps from the tumor" ,
        is_optional => 0,
        is_input => 1},

    normal_snp_file => {
        is => 'Text',
        doc => "all snps from the normal" ,
        is_optional => 0,
        is_input => 1},

    file_source => {
        is => 'Text',
        doc => "source of the bam files",
        is_optional => 1,
        is_input => 1,
        default =>"dbGap" },

    build_dir => {
        is => 'Text',
        doc => "Build directory",
        is_optional => 0,
        is_input => 1},

    dbsnp_file => {
        is => 'Text',
        doc => "dbsnp File - if specified, will label dbSNP sites" ,
        is_optional => 1,
        is_input => 1,
        default => ""},

    individual_id => {
        is => 'Text',
        doc => "Individual ID",
        is_optional => 0,
        is_input => 1},


    center => {
        is => 'Text',
        doc => "Genome center name (WUSTL, Broad, Baylor)" ,
        is_optional => 1,
        default => "WUSTL",
        is_input => 1},


    chrom => {
        is => 'Text',
        doc => "do only this chromosome" ,
        is_optional => 1,
        default => "",
        is_input => 1},

    skip_header => {
        is => 'Boolean',
        is_optional => 1,
        is_input => 1,
        default => 0,
        doc => 'enable this to skip header output - useful for doing individual chromosomes. Note that the output will be appended to the output file if this is enabled.',
    },

    filterOut => {
        is => 'Text',
        doc => "file containing SNVs to label as filtered out in format: filterName:description:file,filterName2,description2,file2,..." ,
        is_optional => 1,
        default => "",
        is_input => 1},


    filterPass => {
        is => 'Text',
        doc => "file containing SNVs that pass filters (all SNVs not contained in the file will be marked filtered out) in format:  filterName:description:file,filterName2,description2:file2,... " ,
        is_optional => 1,
        default => "",
        is_input => 1},


    genome_build => {
        is => 'Text',
        doc => "Reference genome build - currently only supports b36" ,        
        is_optional => 0,
        is_input => 1},

    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Generate Vcf File"
}


sub help_synopsis {
    <<'HELP';
Generate a VCF File
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Parses the relevant files and creates a VCF containing all the SNVs. This includes those that fail filters (noted in the FILTER field).
HELP
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    my $output_file = $self->output_file;
    my $tumor_bam = $self->tumor_bam_file;
    my $normal_bam = $self->normal_bam_file;
    my $file_source = $self->file_source;
    my $build_dir = $self->build_dir;
    my $individual_id = $self->individual_id;
    my $center = $self->center;
    my $genome_build = $self->genome_build;
    my $dbsnp_file = $self->dbsnp_file;
    my $chrom = $self->chrom;
    my $tumor_snp_file = $self->tumor_snp_file;
    my $normal_snp_file = $self->normal_snp_file;
    my $skip_header = $self->skip_header;
    my $filterPass = $self->filterPass;
    my $filterOut = $self->filterOut;

    my $analysis_profile = "samtools pileup and/or somatic-sniper";

###########################################################################

#############################################################################
    sub print_header{
        my ($tumor_bam, $normal_bam, $center, $genome_build, $individual_id, $file_source, $analysis_profile, $output_file) = @_;

        open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
        my $reference;
        my $seqCenter;
        my $file_date = localtime();


        #fix this to support build 37 when necessary
        if ($genome_build ne "36"){
            die("reference paths need to be added for other builds before using")
        }


        #center-specific lines:
        if ($center eq "WUSTL"){
            $seqCenter = "genome.wustl.edu";
            $reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36_WUGSC_variant.fa.gz";
        }
        elsif($center eq "Broad"){
            $seqCenter = "broad.mit.edu";
            $reference="ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36-HG18_Broad_variant.fa.gz";
        }
        elsif($center eq "Baylor"){
            $seqCenter = "bcm.edu";
            $reference="ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36_BCCAGSC_variant.fa.gz";
        }


        #which type of sequencing was it?
        my $seqType = "";
        $seqType = "-whole" if($tumor_bam =~ /whole/i);
        $seqType = "-solid" if($tumor_bam =~ /solid/i);
        $seqType = "-illumina" if($tumor_bam =~ /illumina/i);

        print OUTFILE "##fileformat=VCFv4.1" . "\n";
#	print OUTFILE "##fileDate=" . strftime("%Y/%m/%d %H:%M:%S\n", localtime) . "\n";
        print OUTFILE "##fileDate=" . $file_date . "\n";

        print OUTFILE "##reference=$reference" . "\n";
        print OUTFILE "##phasing=none" . "\n";
        print OUTFILE "##INDIVIDUAL=$individual_id" . "\n";

        #first normal
        print OUTFILE "##SAMPLE=<ID=NORMAL,file=" . $normal_bam . ",SeqCenter=" . $seqCenter . ",Type=normal_DNA>" . "\n";

        #then tumor
        print OUTFILE "##SAMPLE=<ID=TUMOR,file=" . $tumor_bam . ",SeqCenter=" . $seqCenter . ",Type=tumor_DNA>" . "\n";

        # info lines
        print OUTFILE "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 130\">" . "\n";
#	print OUTFILE "##INFO=<ID=VT,Number=1,Type=String,Description=\"Somatic variant type\">" . "\n";

        #all the filter info
        print OUTFILE "##FILTER=<ID=PASS,Description=\"Passed all filters\">" . "\n";
        print OUTFILE "##FILTER=<ID=snpfilter,Description=\"Failed Maq inspired SNPFilter - Discard\">" . "\n";
        print OUTFILE "##FILTER=<ID=loh,Description=\"Loss of Heterozygosity filter - Discard\">" . "\n";
        print OUTFILE "##FILTER=<ID=ceuyri,Description=\"Novel event filter (CEU and YRI) - Discard\">" . "\n";
        print OUTFILE "##FILTER=<ID=dbsnp,Description=\"found in dbSNP - Discard\">" . "\n";
        print OUTFILE "##FILTER=<ID=sniper,Description=\"not called as potentially somatic by somatic sniper - Potentially germline\">" . "\n";

        #format info
        print OUTFILE "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" . "\n";
        print OUTFILE "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" . "\n";
        print OUTFILE "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">" . "\n";
        #print OUTFILE "##FORMAT=<ID=BQ,Number=1,Type=Integer,Description=\"Average Base Quality corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
        #print OUTFILE "##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Average Mapping Quality corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
        #print OUTFILE "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
        print OUTFILE "##FORMAT=<ID=VAS,Number=1,Type=Integer,Description=\"Variant  Status relative to non-adjacent normal 0=Wildtype, 1=Germline, 2=Somatic, 3=LOH, 4=Post_Transcriptional_Modification, 5=Undefined\">" . "\n";
        print OUTFILE "##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description=\"Quality score - SomaticSniper score\">" . "\n";

        #column header:
        print OUTFILE  "#" . join("\t", ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR")) . "\n";
        OUTFILE->close();
    }





###################################################################

    #everything is hashed by chr:position, with subhashes corresponding to
    #tumor and normal samples, then the various VCF fields
    my %allSnvs;

    #First, read in the unfiltered files to get all of the calls
    #tumor calls from samtools
    my $inFh = IO::File->new( "$tumor_snp_file" ) || die "can't open file\n";

    while( my $line = $inFh->getline ){

        chomp($line);
        my @col = split("\t",$line);

        #if we do this on a per-chrom process (for huge files)
        unless ($chrom eq ""){
            next if($col[0] ne $chrom)
        }

        my $chr = $col[0];
        #replace X and Y for sorting
        $chr = "23" if $col[0] eq "X";
        $chr = "24" if $col[0] eq "Y";
        $chr = "25" if $col[0] eq "MT";
        my $id = $chr . ":" . $col[1];

        #skip MT and NT chrs
        #next if $col[0] =~ /^MT/;
        next if $col[0] =~ /^NT/;

        $allSnvs{$id}{"chrom"} = $col[0];
        $allSnvs{$id}{"pos"} = $col[1];

        #replace ambiguous/IUPAC bases with N in ref
        $col[2] =~ s/[^ACGTN]/N/g;

        #get all the alleles together (necessary for the GT field)
        my @allAlleles = ($col[2]);
        my @varAlleles;
        my @tmp = Genome::Info::IUB->iub_to_bases($col[3]);

        #only add non-reference alleles to the alt field
        foreach my $alt (@tmp){
            unless ($alt eq $col[2]){
                push(@allAlleles,$alt);
                push(@varAlleles,$alt);
            }
        }

        $allSnvs{$id}{"ref"} = $col[2];
        $allSnvs{$id}{"alt"} = join(",",@varAlleles);

        #add the ref and alt alleles' positions in the allele array to the GT field
        $allSnvs{$id}{"normal"}{"GT"} = ".";
        $allSnvs{$id}{"tumor"}{"GT"} = genGT($col[3],@allAlleles);

        #genotype quality (consensus quality)
        $allSnvs{$id}{"normal"}{"GQ"} = ".";
        $allSnvs{$id}{"tumor"}{"GQ"} = $col[4];

        #total read depth
        $allSnvs{$id}{"normal"}{"DP"} = ".";
        $allSnvs{$id}{"tumor"}{"DP"} = $col[7];

        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"MQ"} = ".";
        $allSnvs{$id}{"tumor"}{"MQ"} = ".";

        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"BQ"} = ".";
        $allSnvs{$id}{"tumor"}{"BQ"} = ".";

        # #allele depth
        $allSnvs{$id}{"normal"}{"AD"} =  ".";
        $allSnvs{$id}{"tumor"}{"AD"} =  ".";

        # vas -- Variant  Status relative to non-adjacent normal
        # 0=Wildtype, 1=Germline, 2=Somatic, 3=LOH, 4=Post_Transcriptional_Modification, 5=Undefined
        $allSnvs{$id}{"normal"}{"VAS"} = ".";
        $allSnvs{$id}{"tumor"}{"VAS"} = ".";

        #vaq
        $allSnvs{$id}{"normal"}{"VAQ"} = ".";
        $allSnvs{$id}{"tumor"}{"VAQ"} = ".";

        # #vls
        # $allSnvs{$id}{"normal"}{"VLS"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLS"} = ".";

        # #vlq
        # $allSnvs{$id}{"normal"}{"VLQ"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLQ"} = ".";

        $allSnvs{$id}{"filter"} = "sniper";
    }
    $inFh->close();



#-------------------------------------------
    #now normal calls from samtools
    my $inFh3 = IO::File->new( "$normal_snp_file" ) || die "can't open file\n";
    while( my $line = $inFh3->getline ){

        chomp($line);
        my @col = split("\t",$line);

        #if we do this on a per-chrom process (for huge files)
        unless ($chrom eq ""){
            next if($col[0] ne $chrom)
        }

        my $chr = $col[0];
        #replace X and Y for sorting
        $chr = "23" if $col[0] eq "X";
        $chr = "24" if $col[0] eq "Y";
        $chr = "25" if $col[0] eq "MT";
        my $id = $chr . ":" . $col[1];

        #skip MT and NT chrs
        #next if $col[0] =~ /^MT/;
        next if $col[0] =~ /^NT/;

        #replace ambiguous/IUPAC bases with N in ref
        $col[2] =~ s/[^ACGTN]/N/g;
        
        #if we didn't also see this pos as a tumor snv
        if (!(exists($allSnvs{$id}))){
            $allSnvs{$id}{"chrom"} = $col[0];
            $allSnvs{$id}{"pos"} = $col[1];
        }

        my @allAlleles = ($col[2]);
        my @varAlleles;
        # if exists, we have to consider both tumor and variant alleles for gt positions
        #retrieve the ref/alt from the tumor call
        if (exists($allSnvs{$id}{"alt"})){
            my @tmp = split(",", $allSnvs{$id}{"alt"});
            @allAlleles = (@allAlleles, @tmp);
            @varAlleles = @tmp;
        }

        my @tmp = Genome::Info::IUB->iub_to_bases($col[3]);
        #only add non-reference alleles to the alt field
        foreach my $alt (@tmp){
            unless (grep $_ eq $alt, @allAlleles){
                push(@allAlleles,$alt);
                push(@varAlleles,$alt);
            }
        }

        $allSnvs{$id}{"ref"} = $col[2];
        #replace alt in case we have added more alleles
        $allSnvs{$id}{"alt"} = join(",",@varAlleles);

        #add the ref and alt alleles' positions in the allele array to the GT field
        $allSnvs{$id}{"normal"}{"GT"} = genGT($col[3],@allAlleles);

        #genotype quality (consensus quality)
        $allSnvs{$id}{"normal"}{"GQ"} = $col[4];
        #total read depth
        $allSnvs{$id}{"normal"}{"DP"} = $col[7];
        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"MQ"} = ".";
        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"BQ"} = ".";

        #dot out same three in tumor if we don't have vals there already
        if(!(exists($allSnvs{$id}))){
            $allSnvs{$id}{"tumor"}{"GT"} = ".";
            $allSnvs{$id}{"tumor"}{"GQ"} = ".";
            $allSnvs{$id}{"tumor"}{"DP"} = ".";
            $allSnvs{$id}{"tumor"}{"MQ"} = ".";
            $allSnvs{$id}{"tumor"}{"BQ"} = ".";
        }

        #allele depth
        $allSnvs{$id}{"normal"}{"AD"} =  ".";
        $allSnvs{$id}{"tumor"}{"AD"} =  ".";

        # vas -- Variant  Status relative to non-adjacent normal
        # 0=Wildtype, 1=Germline, 2=Somatic, 3=LOH, 4=Post_Transcriptional_Modification, 5=Undefined
        $allSnvs{$id}{"normal"}{"VAS"} = ".";
        $allSnvs{$id}{"tumor"}{"VAS"} = ".";

        #vaq
        $allSnvs{$id}{"normal"}{"VAQ"} = ".";
        $allSnvs{$id}{"tumor"}{"VAQ"} = ".";

        # #vls
        # $allSnvs{$id}{"normal"}{"VLS"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLS"} = ".";
        # #vlq
        # $allSnvs{$id}{"normal"}{"VLQ"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLQ"} = ".";

        $allSnvs{$id}{"filter"} = "sniper";
    }
    $inFh->close();



#-------------------------------------------

    #now the sniper calls
    $inFh = IO::File->new( "$build_dir/snv_sniper_snp.csv" ) || die "can't open file\n";
    while( my $line = $inFh->getline ){

        chomp($line);
        my @col = split("\t",$line);

        #if we do this on a per-chrom process (for huge files)
        unless ($chrom eq ""){
            next if($col[0] ne $chrom)
        }

        my $chr = $col[0];
        #replace X and Y for sorting
        $chr = "23" if $col[0] eq "X";
        $chr = "24" if $col[0] eq "Y";
        $chr = "25" if $col[0] eq "MT";
        my $id = $chr . ":" . $col[1];

        #skip MT and NT chrs
        #next if $col[0] =~ /^MT/;
        next if $col[0] =~ /^NT/;

        if (!(exists($allSnvs{$id}))){
            #add whole new item to hash	    
            $allSnvs{$id}{"chrom"} = $col[0];
            $allSnvs{$id}{"pos"} = $col[1];
        }    

        #replace ambiguous/IUPAC bases with N in ref
        $col[2] =~ s/[^ACGTN]/N/g;

        #just replace anything from samtools with the presumably better sniper call.
        my @allAlleles = $col[2];
        my @varAlleles;
        my @tmp = Genome::Info::IUB->iub_to_bases($col[3]);
        #only add non-reference alleles to the alt field
        foreach my $alt (@tmp){
            unless ($alt eq $col[2]){
                push(@allAlleles,$alt);
                push(@varAlleles,$alt);
            }
        }

        $allSnvs{$id}{"ref"} = $col[2];
        $allSnvs{$id}{"alt"} = join(",",@varAlleles);

        #add the ref and alt alleles' positions in the allele array to the GT field
        $allSnvs{$id}{"normal"}{"GT"} = genGT($col[2],@allAlleles);
        $allSnvs{$id}{"tumor"}{"GT"} = genGT($col[3],@allAlleles);


        #genotype quality
        $allSnvs{$id}{"tumor"}{"GQ"} = $col[6];
        $allSnvs{$id}{"normal"}{"GQ"} = ".";

        #total read depth
        $allSnvs{$id}{"normal"}{"DP"} = $col[9];
        $allSnvs{$id}{"tumor"}{"DP"} = $col[8];

        #avg base quality ref/var
        $allSnvs{$id}{"tumor"}{"BQ"} =  ".";
        $allSnvs{$id}{"normal"}{"BQ"} =  ".";

        #avg mapping quality ref/var
        $allSnvs{$id}{"tumor"}{"MQ"} =   ".";
        $allSnvs{$id}{"normal"}{"MQ"} =  ".";

        #allele depth
        $allSnvs{$id}{"normal"}{"AD"} =  ".";
        $allSnvs{$id}{"tumor"}{"AD"} =  ".";	

        #vas
        $allSnvs{$id}{"normal"}{"VAS"} = 0;
        $allSnvs{$id}{"tumor"}{"VAS"} = 2;

        #vaq
        $allSnvs{$id}{"normal"}{"VAQ"} = ".";
        $allSnvs{$id}{"tumor"}{"VAQ"} = $col[4];

        # #vls
        # $allSnvs{$id}{"normal"}{"VLS"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLS"} = ".";

        # #vlq
        # $allSnvs{$id}{"normal"}{"VLQ"} = ".";
        # $allSnvs{$id}{"tumor"}{"VLQ"} = ".";

        $allSnvs{$id}{"filter"} = "";
    }
    $inFh->close();


#-------------------------------------------

#Next, go through all the filtered files, match up the snps,
#and add a label to the filter field if it's removed

    sub addFilterInfo{
        my ($filename,$filtername,$snvHashRef,$build_dir) = @_;

        #read in all the sites that passed the filter
        my %passingSNVs;

        my $inFh2 = IO::File->new( "$build_dir/$filename" ) || die "can't open file - $build_dir/$filename\n";
        while( my $line = $inFh2->getline )
        {
            chomp($line);
            my @col = split("\t",$line);
            my $id = $col[0] . ":" . $col[1];

            $passingSNVs{$id} = 1;
        }
        $inFh2->close();

        #check each stored SNV
        foreach my $key (keys( %{$snvHashRef} )){
            #if it did not pass this filter
            unless(exists($passingSNVs{$key})){
                #and hasn't already been filtered out
                if(!(exists($snvHashRef->{$key}{"filter"})) || ($snvHashRef->{$key}{"filter"} eq "")){ 
                    #add the filter name
                    $snvHashRef->{$key}{"filter"} = $filtername;
                }
            }
        }
    }


#---------------------------------------------
##ADD FILTER INFO HERE
##like this: addFilterInfo(filename filtername, sniperSNPhash)
    addFilterInfo("sfo_snp_filtered.csv","snpfilter",\%allSnvs, $build_dir);
    addFilterInfo("noloh.csv","loh",\%allSnvs, $build_dir);
    addFilterInfo("ceu_yri_filtered.csv","ceuyri",\%allSnvs, $build_dir);
    addFilterInfo("dbsnp_filtered.csv","dbsnp",\%allSnvs, $build_dir);


    # sub dedupFilterNames{
    # 	my ($names1,$names2) = @_;
    # 	my @n1 = split(",",$names1);
    # 	my @n2 = split(",",$names2);
    # 	return(join(";",uniq(sort(@n1,@n2))))
    # }


#---------------------------------------------
    ## add DBsnp labels, if --dbsnp is specified
    if ($dbsnp_file ne ""){

        print STDERR "adding dbSNP info - this will take a few minutes\n";
        my $inFh = IO::File->new( $dbsnp_file ) || die "can't open file\n";
        while( my $line = $inFh->getline )
        {
            unless($line =~ /^#/){
                chomp($line);
                my @fields = split("\t",$line);

                $fields[1] =~ s/chr//;

                #replace X and Y for sorting
                my $chr = $fields[1];
                $chr = "23" if $chr eq "X";
                $chr = "24" if $chr eq "Y";

                #ucsc is zero-based, so we adjust
                my $key = $chr . ":" . ($fields[2]+1);
                #if the line matches this dbsnp position
                if(exists($allSnvs{$key})){
                    #and the alleles match
                    if (($allSnvs{$key}{"alt"} . "/" . $allSnvs{$key}{"ref"} eq $fields[9])){
                        #note the match in the info field
                        if(exists($allSnvs{$key}{"info"})){
                            $allSnvs{$key}{"info"} = $allSnvs{$key}{"info"} . ";";
                        } else {
                            $allSnvs{$key}{"info"} = "";
                        }
                        $allSnvs{$key}{"info"} = $allSnvs{$key}{"info"} . "DB";

                        #add to id field
                        if(exists($allSnvs{$key}{"id"})){
                            $allSnvs{$key}{"id"} = $allSnvs{$key}{"id"} . ";";
                        } else {
                            $allSnvs{$key}{"id"} = "";
                        }
                        $allSnvs{$key}{"id"} = $allSnvs{$key}{"id"} . $fields[4];


#			#if the filter shows a pass, remove it and add dbsnp
#			if($allSnvs{$key}->{FILTER} eq "PASS"){
#			    $allSnvs{$key}->{FILTER} = "dbSNP";
#			} else { #add dbsnp to the list
#			    $allSnvs{$key}->{FILTER} = $allSnvs{$key}->{FILTER} . ",dbSNP";
#			}

                    }
                }
            }
        }
    }


#---------------------------------------------
    sub print_body{
        my ($output_file,$snvHash) = @_;

        open(OUTFILE, ">>$output_file") or die "Can't open output file: $!\n";
        my %snvhash = %{$snvHash};

        #sort by chr, start for clean output
        sub keySort{
            my($x,$y) = @_;
            my @x1 = split(":",$x);
            my @y1 = split(":",$y);
            return($x1[0] <=> $y1[0] || $x1[1] <=> $y1[1])
        }
        my @sortedKeys = sort { keySort($a,$b) } keys %snvhash;

        foreach my $key (@sortedKeys){
            my @outline;
            push(@outline, $snvhash{$key}{"chrom"});
            push(@outline, $snvhash{$key}{"pos"});


            #ID
            if (exists($snvhash{$key}{"id"})){
                push(@outline, $snvhash{$key}{"id"});
            } else {
                push(@outline, ".");
            }

            #ref/alt
            push(@outline, $snvhash{$key}{"ref"});
            push(@outline, $snvhash{$key}{"alt"});

            #QUAL
            if (exists($snvhash{$key}{"qual"})){
                push(@outline, $snvhash{$key}{"qual"});
            } else {
                push(@outline, ".");
            }

            #FILTER
            if (exists($snvhash{$key}{"filter"}) && $snvhash{$key}{"filter"} ne ""){
                push(@outline, $snvhash{$key}{"filter"});
            } else {
                push(@outline, "PASS");
            }

            #INFO
            push(@outline, ".");

            #FORMAT
            # push(@outline, "GT:GQ:DP:BQ:MQ:AD:VAS:VAQ:VLS:VLQ"); 
            #push(@outline, "GT:GQ:DP:BQ:MQ:AD:VAS:VAQ");
            push(@outline, "GT:GQ:DP:VAS:VAQ");

            my @normalFormat;
            my @tumorFormat;

            #my @fields = ("GT","GQ","DP","BQ","MQ","AD","VAS","VAQ");
            my @fields = ("GT","GQ","DP","VAS","VAQ");
            #collect format fields
            foreach my $field (@fields){
                if(exists($snvhash{$key}{"normal"}{$field})){
                    push(@normalFormat, $snvhash{$key}{"normal"}{$field});
                } else {
                    push(@normalFormat,".")
                }
                if(exists($snvhash{$key}{"tumor"}{$field})){
                    push(@tumorFormat, $snvhash{$key}{"tumor"}{$field});
                } else {
                    push(@tumorFormat,".")
                }

            }
            push(@outline, join(":",@normalFormat));
            push(@outline, join(":",@tumorFormat));

            print OUTFILE join("\t",@outline) . "\n";
        }
    }

#----------------------------------
    unless ($skip_header){
        print_header($tumor_bam, $normal_bam, $center, $genome_build, $individual_id, $file_source, $analysis_profile,$output_file);
    }
    print_body($output_file, \%allSnvs);
    return 1;
}
