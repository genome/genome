package Genome::Model::Tools::Vcf::VcfMakerSamtools;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);


class Genome::Model::Tools::Vcf::VcfMakerSamtools {
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
            doc => "Reference genome build" ,
            is_optional => 1,
            default => "36",
        },
        
        samtools_file => {
            is => 'Text',
            doc => "samtools output file" ,
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

        dbsnp_file => {
            is => 'Text',
            doc => "dbsnp File - if specified, will label dbSNP sites",
            is_optional => 1,
            is_input => 1,
            default => "",
        },

        ],
};


sub help_brief {                            # keep this to just a few words <---
    "Generate Vcf File from samtools output"
}


sub help_synopsis {
    <<'HELP';
    Generate a VCF File from samtools output
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
    my $genome_build = $self->genome_build;
    my $chrom = $self->chrom;
    my $skip_header = $self->skip_header;
    my $samtools_file = $self->samtools_file;
    my $sample_id = $self->sample_id;
    my $type = $self->type;
    my $dbsnp_file = $self->dbsnp_file;


    if(($type ne "snv") && ($type ne "indel")){
        die("\"type\" parameter must be one of \"snv\" or \"indel\"");
    }

    my $genome_build_fasta;
    if ($genome_build =~ m/36/) {
        my $genome_build_fasta_object= Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $genome_build_fasta = $genome_build_fasta_object->data_directory . "/all_sequences.fa";
    }
    elsif ($genome_build =~ m/37/) {
        my $genome_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $genome_build_fasta = $genome_build_fasta_object->data_directory . "/all_sequences.fa";
    }
    else {
        die "Please specify either build 36 or 37";
    }

###########################################################################
# subs

    sub convertIub{
        my ($base) = @_;
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
    };


    #------------------------
    sub genGT{
        my ($base, @alleles) = @_;
        my @bases = split(",",convertIub($base));
        if (@bases > 1){
            my @pos;
            push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
            push(@pos, (firstidx{ $_ eq $bases[1] } @alleles));
            return(join("/", sort(@pos)));
        } else { #only one base
            my @pos;
            push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
            push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
            return(join("/", sort(@pos)));
        }
    }

    #-------------------------
    #get preceding base using samtools faidx
    sub getPrecedingBase{
        my ($chr,$pos,$genome_build_fasta) = @_;
        my $base = `samtools faidx $genome_build_fasta $chr:$pos-$pos | tail -n 1`;
        chomp($base);
        return($base)
    }


    #-------------------------
    sub print_header{
        my ($genome_build, $sample_id, $output_file) = @_;

        open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
        my $seqCenter;
        my $file_date = localtime();
        my $reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36_BCCAGSC_variant.fa.gz";

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
        print OUTFILE  "#" . join("\t", ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","$sample_id")) . "\n";
        OUTFILE->close();
    }



#-------------------------------------------
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
            if (exists($snvhash{$key}{"INFO"})){
                push(@outline, $snvhash{$key}{"INFO"});
            } else {
                push(@outline, ".");
            }

            #FORMAT
            push(@outline, "GT:GQ:DP:BQ:MQ:AD:FA:VAQ");

            my @format;

            my @fields = ("GT","GQ","DP","BQ","MQ","AD","FA","VAQ");
            #collect format fields
            foreach my $field (@fields){
                if(exists($snvhash{$key}{$field})){
                    push(@format, $snvhash{$key}{$field});
                } else {
                    push(@format,".")
                }

            }
            push(@outline, join(":",@format));

            print OUTFILE join("\t",@outline) . "\n";
        }
    }


###################################################################
# actually do the parsing here
    sub samtoolsRead{
        my ($samtools_file, $chrom, $type, $genome_build_fasta) = @_;

        #everything is hashed by chr:position, with subhashes corresponding to
        #samples, then the various VCF fields
        my %allSnvs;

        #First, read in the unfiltered files to get all of the calls
        #tumor calls from samtools
        my $inFh = IO::File->new( "$samtools_file" ) || die "can't open file\n";

        while( my $line = $inFh->getline ){

            chomp($line);
            my @col = split("\t",$line);

            #if we do this on a per-chrom process (for huge files)
            unless ($chrom eq ""){
                next if($col[0] ne $chrom)
            }

            my $chr = $col[0];
            #replace X and Y for sorting
            $ chr = "23" if $col[0] eq "X";
            $chr = "24" if $col[0] eq "Y";
            $chr = "25" if $col[0] eq "MT";

            #replace ambiguous/IUPAC bases with N in ref
            $col[3] =~ s/[^ACGTN\-]/N/g;
            
            #get slash-sep ids for variant
            my $altvar = join("/", split(",",convertIub($col[4])));
            my $id = $chr . ":" . $col[1] . ":" . $col[2] . ":" . $col[3] . ":" . $altvar;

            #skip MT and NT chrs
            #next if $col[0] =~ /^MT/;
            next if $col[0] =~ /^NT/;

            $allSnvs{$id}{"chrom"} = $col[0];
            $allSnvs{$id}{"pos"} = $col[1];


            #handle snv genotype calls
            if ($type eq "snv"){        
                #get all the alleles together (necessary for the GT field)
                my @allAlleles = $col[3];
                my @varAlleles;
                my @tmp = split(",",convertIub($col[4]));
                
                #only add non-reference alleles to the alt field
                foreach my $alt (@tmp){
                    unless ($alt eq $col[3]){
                        push(@allAlleles,$alt);
                        push(@varAlleles,$alt);
                    }
                }

                $allSnvs{$id}{"ref"} = $col[3];
                $allSnvs{$id}{"alt"} = join(",",@varAlleles);

                #add the ref and alt alleles' positions in the allele array to the GT field
                $allSnvs{$id}{"GT"} = genGT($col[4],@allAlleles);
                $allSnvs{$id}{"INFO"} = "VT=SNP";




                #handle indel genotype calls
            } elsif ($type eq "indel"){
                #add the preceding base as an anchor position
                my $pbase = getPrecedingBase($col[0],$col[1],$genome_build_fasta);



                #insertion
                if ($col[3] eq "-"){  
                    $allSnvs{$id}{"ref"} = $pbase;
                    $allSnvs{$id}{"alt"} = $pbase . $col[4];
                    
                    #deletion
                } elsif ($col[4] eq "-"){
                    $allSnvs{$id}{"ref"} = $pbase . $col[3];
                    $allSnvs{$id}{"alt"} = $pbase;
                    #confusion
                } else {
                    die("this isn't an insertion or deletion - what is it?\n$line");
                }

                $allSnvs{$id}{"GT"} = "0/1";

                #we have a second indel at that spot, with no associated
                #quality information - just print an error for now
                if (($col[9] ne "*") && ($col[10] ne "*")){
                    print STDERR "second indel at this pos with no qual information is being skipped:\n$line\n";
                }
                $allSnvs{$id}{"INFO"} = "VT=INDEL";

            }


            #genotype quality (consensus quality)
            $allSnvs{$id}{"GQ"} = $col[5];

            #total read depth
            $allSnvs{$id}{"DP"} = $col[8];

            #avg mapping quality ref/var
            $allSnvs{$id}{"MQ"} = ".";

            #avg mapping quality ref/var
            $allSnvs{$id}{"BQ"} = ".";

            # #allele depth
            $allSnvs{$id}{"AD"} =  ".";
            
            #fraction of reads supporting alt
            $allSnvs{$id}{"FA"} =  ".";

            #vaq
            $allSnvs{$id}{"VAQ"} = ".";

        }
        $inFh->close();

        return %allSnvs;
    }

#---------------------------------------------
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
    

#----------------------------------
    my %samtools_hash = samtoolsRead($samtools_file, $chrom, $type, $genome_build_fasta);

    
    unless ($skip_header){
        print_header($genome_build, $sample_id, $output_file);
    }

    ## add DBsnp labels, if --dbsnp is specified
    if ($dbsnp_file ne ""){
        addDbSnp($dbsnp_file, $chrom, \%samtools_hash)
    }

    print_body($output_file, \%samtools_hash);

    return 1;
}
