package Genome::Model::Tools::Vcf::VcfMakerSomaticBase;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Vcf::VcfMakerSomaticBase {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "List of mutations in Vcf format",
        },
        chrom => {
            is => 'Text',
            doc => "do only this chromosome" ,
            is_optional => 1,
            default => "",
        },
        skip_header => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this to skip header output - useful for doing individual chromosomes. Note that the output will be appended to the output file if this is enabled.',
        },
        genome_build => {
            is => 'Text',
        },
        input_file => {
            is => 'Text',
            doc => "variant detector output file" ,
            is_optional => 0,
            is_input => 1,
        },
        type => {
            is => 'Text',
            doc => "type of variant calls - one of \"snv\" or \"indel\"" ,
            is_optional => 0,
            is_input => 1,
        },
        sample_id => {
            is => 'Text',
            doc => "unique sample id",
            is_optional => 0,
            is_input => 1,
        },
        standard_chroms => {
            default=>0,
            doc=> "set to 1 if you only want 1..2,X,Y,MT"
        },
        dbsnp_file => {
            is => 'Text',
            doc => "dbsnp File - if specified, will label dbSNP sites",
            is_optional => 1,
            is_input => 1,
            default => "",
        },
        seq_center => {
            is => 'Text',
            doc => "Center that did the sequencing (WUSTL or BROAD)" ,
            is_optional => 1,
            default => "WUSTL",
        },
        cp_score_to_qual => {
            is => 'Boolean',
            doc => "copy the somatic score to the qual field for Mutation WG comparisons" ,
            is_optional => 1,
            default => 0,
            is_input => 1
        },
    ],
};

sub print_header{
    my ($genome_build, $sample_id, $output_file, $seq_center) = @_;

    open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
    my $file_date = localtime();

    unless($genome_build eq "36"){
        die("reference paths only provided for b36 for broad and wustl. Update the tool with the appropriate paths for other builds");
    }

    my $reference;
    if ($seq_center eq "WUSTL"){
        $reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36_BCCAGSC_variant.fa.gz";
    } elsif ($seq_center eq "BROAD"){
        $reference="ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36-HG18_Broad_variant.fa.gz";
    } else {
        die ("only have references hardcoded for WUSTL and BROAD");
    }

    print OUTFILE "##fileformat=VCFv4.0" . "\n";
    print OUTFILE "##fileDate=" . $file_date . "\n";
    print OUTFILE "##reference=$reference" . "\n";
    print OUTFILE "##phasing=none" . "\n";
    print OUTFILE "##SAMPLE=$sample_id" . "\n";

#format info
    print OUTFILE "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" . "\n";
    print OUTFILE "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" . "\n";
    print OUTFILE "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">" . "\n";
    print OUTFILE "##FORMAT=<ID=BQ,Number=1,Type=Integer,Description=\"Average Base Quality corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
    print OUTFILE "##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Average Mapping Quality corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
    print OUTFILE "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth corresponding to alleles 0/1/2/3... after software and quality filtering\">" . "\n";
    print OUTFILE "##FORMAT=<ID=FA,Number=1,Type=Float,Description=\"Fraction of reads supporting ALT\">" . "\n";
    print OUTFILE "##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description=\"Variant Quality\">" . "\n";

#INFO
    print OUTFILE "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type\">" . "\n";

#column header:
    print OUTFILE  "#" . join("\t", ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","PRIMARY")) . "\n";
    OUTFILE->close();
}

sub print_body{
    my ($output_file,$snvHash,$cp_score_to_qual) = @_;

    open(OUTFILE, ">>$output_file") or die "Can't open output file: $!\n";
    my %snvhash = %{$snvHash};

    #sort by chr, start for clean output
    sub keySort{
        my($x,$y) = @_;
        my @x1 = split(":",$x);
        my @y1 = split(":",$y);
        return($x1[0] cmp $y1[0] || $x1[1] <=> $y1[1])
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
            if ($cp_score_to_qual){
                push(@outline,$snvhash{$key}{"tumor"}{"VAQ"});
            } else {
                push(@outline, ".");
            }
        }

        #FILTER
        if (exists($snvhash{$key}{"filter"}) && $snvhash{$key}{"filter"} ne ""){
            push(@outline, $snvhash{$key}{"filter"});
        } else {
            push(@outline, "PASS");
        }

        #INFO
        if (exists($snvhash{$key}{"INFO"})){
            push(@outline, $snvhash{$key}{"INFO"});
        } else {
            push(@outline, ".");
        }

        #FORMAT
        push(@outline, "GT:GQ:DP:BQ:MQ:AD:FA:VAQ");

        my @fields = ("GT","GQ","DP","BQ","MQ","AD","FA","VAQ");

        #collect format fields
        foreach my $type ("normal","tumor"){
            my @format;
            foreach my $field (@fields){
                if(exists($snvhash{$key}{$type}{$field})){
                    push(@format, $snvhash{$key}{$type}{$field});
                } else {
                    push(@format,".")
                }

            }
            push(@outline, join(":",@format));
        }

        print OUTFILE join("\t",@outline) . "\n";
    }
}

sub runGetPrecedingBase {
    my ($chr, $pos, $cmd) = @_;
    my $base = `$cmd`;
    chomp($base);
    return($base);
}

sub readVariantFile{
    my ($file, $chrom, $type, $standard_chroms) = @_;

    #everything is hashed by chr:position, with subhashes corresponding to
    #samples, then the various VCF fields
    my %allSnvs;

    #First, read in the unfiltered files to get all of the calls
    #tumor calls from sniper
    my $inFh = IO::File->new( "$file" ) || die "can't open file\n";

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

        #replace ambiguous/IUPAC bases with N in ref
        $col[2] =~ s/[^ACGTN\-]/N/g;

        my $id = $chr . ":" . $col[1] . ":" . $col[2] . ":" . $col[3];

        #skip non-normal chrs
        #next if $col[0] =~ /^NT/;

        next if ($standard_chroms && !($col[0] =~/^[1]?([0-9]|^2[012]|X|Y|MT)$/));

        $allSnvs{$id}{"chrom"} = $col[0];
        $allSnvs{$id}{"pos"} = $col[1];


        #handle snv genotype calls
        if ($type eq "snv"){
            #get all the alleles together (necessary for the GT field)
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
            $allSnvs{$id}{"normal"}{"GT"} = ".";
            $allSnvs{$id}{"tumor"}{"GT"} = genGT($col[3],@allAlleles);
            $allSnvs{$id}{"INFO"} = "VT=SNP";



            #handle indel genotype calls
        } elsif ($type eq "indel"){
            #add the preceding base as an anchor position
            my $pbase = getPrecedingBase($col[0],$col[1]);


            #insertion
            if ($col[2] eq "-"){
                $allSnvs{$id}{"ref"} = $pbase;
                $allSnvs{$id}{"alt"} = $pbase . $col[3];

                #deletion
            } elsif ($col[3] eq "-"){
                $allSnvs{$id}{"ref"} = $pbase . $col[2];
                $allSnvs{$id}{"alt"} = $pbase;
                #confusion
            } else {
                die("this isn't an insertion or deletion - what is it?\n$line");
            }

##???
            $allSnvs{$id}{"GT"} = "0/1";

            #we have a second indel at that spot, with no associated
            #quality information - just print an error for now
            if (($col[8] ne "*") && ($col[9] ne "*")){
                print STDERR "second indel at this pos with no qual information is being skipped:\n$line\n";
            }
            $allSnvs{$id}{"INFO"} = "VT=INDEL";

        }


        #genotype quality (consensus quality)
        $allSnvs{$id}{"normal"}{"GQ"} = ".";
        $allSnvs{$id}{"tumor"}{"GQ"} = $col[6];

        #total read depth
        $allSnvs{$id}{"normal"}{"DP"} = $col[9];
        $allSnvs{$id}{"tumor"}{"DP"} = $col[8];

        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"MQ"} = ".";
        $allSnvs{$id}{"tumor"}{"MQ"} = $col[6];

        #avg mapping quality ref/var
        $allSnvs{$id}{"normal"}{"BQ"} = ".";
        $allSnvs{$id}{"tumor"}{"BQ"} = $col[7];

        # #allele depth
        $allSnvs{$id}{"normal"}{"AD"} =  ".";
        $allSnvs{$id}{"tumor"}{"AD"} =  ".";

        #fraction of reads supporting alt
        $allSnvs{$id}{"normal"}{"FA"} =  ".";
        $allSnvs{$id}{"tumor"}{"FA"} =  ".";

        #vaq
        $allSnvs{$id}{"normal"}{"VAQ"} = ".";
        $allSnvs{$id}{"tumor"}{"VAQ"} = $col[4];


    }
    $inFh->close();

    return %allSnvs;
}

sub addDbSnp{
    my ($dbsnp_file, $chrom, $allsnvs)= @_;
    my %allSnvs = %{$allsnvs};

    print STDERR "adding dbSNP info - this will take a few minutes\n";
    my $inFh = IO::File->new( $dbsnp_file ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        unless($line =~ /^#/){
            chomp($line);
            my @fields = split("\t",$line);

            $fields[1] =~ s/chr//;

            #skip snps on chrs we're not considering
            if($chrom ne ""){
                next if $fields[1] ne $chrom;
            }

            #replace X and Y for sorting
            my $chr = $fields[1];
            $chr = "23" if $chr eq "X";
            $chr = "24" if $chr eq "Y";

            #ucsc is zero-based, so we adjust
            my $pos = $fields[2]+1;

            my $key = $chr . ":" . $pos . ":" . $fields[7] . ":" . $fields[9];

            #if the line matches this dbsnp position
            if(exists($allSnvs{$key})){
                #and the alleles match
                # #note the match in the info field
                # if(exists($allSnvs{$key}{"info"})){
                #     $allSnvs{$key}{"info"} = $allSnvs{$key}{"info"} . ";";
                # } else {
                #     $allSnvs{$key}{"info"} = "";
                # }
                # $allSnvs{$key}{"info"} = $allSnvs{$key}{"info"} . "DB";

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

1;
