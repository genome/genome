package Genome::Model::Tools::Relationship::MergeAndFixVcfs;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use Genome::Utility::Vcf "open_vcf_file";
use POSIX;
our $DEFAULT_VERSION = '0.02';
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::Relationship::MergeAndFixVcfs {
    is => 'Command',
    has_optional_input => [
    denovo_vcf => {
        is=>'Text',
        is_optional=>1,
        default=>0,
    },
    output_vcf => {
        is=>'Text',
        is_optional=>0,
        is_output=>1,
    },
    standard_vcf => {
        is=>'Text',
        is_optional=>0,
    },
    bgzip => {
        is_optional=>1,
        default=>1,
        doc=>'set this to 0 if you prefer uncompressed',
    },
    ],
    has_param => [
    lsf_resource => {
        is => 'Text',
        default => "-R 'span[hosts=1] rusage[mem=1000] -n 4'",
    },
    lsf_queue => {
        is => 'Text',
        default => 'long',
    },
    ],
};

sub help_brief {
}

sub help_detail {
}

sub execute {
    $DB::single=1;
    my $self=shift;
    my ($denovo_vcf, $standard_vcf) = ($self->denovo_vcf, $self->standard_vcf);
    unless(-s $denovo_vcf) {
        $self->error_message("Denovo VCF $denovo_vcf has no size or is not found. Exiting");
        die;
    }
    unless(-s $standard_vcf) {
        $self->error_message("Standard VCF $standard_vcf has no size or is not found. Exiting");
        die;
    }
    my $fixed_denovo_vcf = $self->fix_vcf_header($denovo_vcf);
    my $fixed_standard_vcf = $self->fix_vcf_header($standard_vcf);
    $DB::single=1;
    my $temp_output_vcf = $self->merge_vcfs($fixed_denovo_vcf, $fixed_standard_vcf);
    $DB::single=1;
    my $output_vcf = $self->output_vcf;
    if($self->bgzip) {
        my $cmd = "bgzip $temp_output_vcf";
        Genome::Sys->shellcmd(cmd=>$cmd); 
        $temp_output_vcf.=".gz";
        `cp $temp_output_vcf $output_vcf`;
    }else { 
        `cp $temp_output_vcf $output_vcf`;
    }
    return 1;
}

sub merge_vcfs {
    my ($self, $denovo_vcf, $standard_vcf) = @_;
    $DB::single=1;
    my ($temp_ofh, $temp_output_filename) = Genome::Sys->create_temp_file();
    $self->write_merged_header($temp_ofh, $denovo_vcf, $standard_vcf);
    my $denovo_sites = $self->objectify_denovo_vcf($denovo_vcf);
    my $standard_fh = Genome::Sys->open_file_for_reading($standard_vcf);
    my ($last_chr, $last_pos) = 0;
    my @denovo_positions;
    my $merged=0;
    while(my $line = $standard_fh->getline) {
        if($line =~m/^#/) { next;}
        if($line =~m/ERROR/) { 
            $self->error_message("Omitting line $line because it contains no useful information\n"); 
            next;
        }
        chomp($line);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
        if($chr ne $last_chr) {
            @denovo_positions = sort {$a <=> $b} keys %{$denovo_sites->{$chr}};
            $last_pos = 0;
        }
        for my $denovo_position (@denovo_positions) {
            if($denovo_position eq $pos) {
                $temp_ofh->print($self->merge_line($line, $denovo_sites->{$chr}->{$pos}));
                delete $denovo_sites->{$chr}->{$pos};
                $merged=1;
            }
            elsif(($denovo_position < $pos) && ($denovo_position > $last_pos)) {
                $temp_ofh->print($self->fix_denovo_vcf_line($denovo_sites->{$chr}->{$denovo_position}));
                delete $denovo_sites->{$chr}->{$denovo_position};
            }
        }
        if($merged) {
            $merged=0;
        }
        else {
            $temp_ofh->print($line . "\n");
        }
        $last_chr = $chr;
        $last_pos = $pos;
    }
    for my $chr (keys %{$denovo_sites}) { #did we fail to use any-- i.e. was standard missing a whole chromosome
        for my $pos (sort {$a <=> $b} keys %{$denovo_sites->{$chr}}) {
            $temp_ofh->print($self->fix_denovo_vcf_line($denovo_sites->{$chr}->{$pos}));
            delete $denovo_sites->{$chr}->{$pos};
        }
    }

    $temp_ofh->close;
    return $temp_output_filename;
}

sub write_merged_header {
    my ($self, $fh, $denovo_vcf, $standard_vcf) = @_;
    my @info_lines_to_add = (qq|##INFO=<ID=DQ,Number=1,Type=Float,Description="De Novo Mutation Quality">|);
    push @info_lines_to_add, qq|##INFO=<ID=DA,Number=.,Type=Integer,Description="De Novo Mutation Allele">|;
    my @format_lines_to_add = (qq|##FORMAT=<ID=DNGL,Number=10,Type=Integer,Description="Denovo Genotype Likelihoods">|);
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGT,Number=1,Type=String,Description="Genotype">|;
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGQ,Number=1,Type=Integer,Description="Genotype Quality">|;
    ####hack to add in header for bingshan's tag that he hasn't added himself
    chomp(my @number_of_AB_tags = `cat $standard_vcf | grep AB`);
    chomp(my @denovo_number_of_AB_tags = `cat $denovo_vcf | grep AB`);
    if(@number_of_AB_tags || @denovo_number_of_AB_tags) {
        my $has_header=0;
        for my $line (@number_of_AB_tags, @denovo_number_of_AB_tags) {
            if($line=~m/^##INFO=<ID=AB,/) {
                $has_header=1;
            }
        }
        unless($has_header) {
            push @info_lines_to_add, qq|##INFO=<ID=AB,Number=1,Type=Float,Description="Allelic Balance">|;
        }
    }
    ######

                   
    my $ifh = Genome::Sys->open_file_for_reading($standard_vcf);
    my ($info_printed, $format_printed)=(0,0);
    while(my $line = $ifh->getline) {
        if($line =~m/^#/) {
            if($line =~m/#INFO/ && !$info_printed) {
                for my $info_line (@info_lines_to_add) {
                    $fh->print($info_line ."\n");
                }
                $info_printed=1;
            }
            if($line =~m/#FORMAT/ && !$format_printed) {
                for my $format_line (@format_lines_to_add) {
                    $fh->print($format_line ."\n");
                }
                $format_printed=1
            }
            $fh->print($line);
        }
        else {
            last;
        }
    }
}



sub merge_line {
    my ($self, $standard_line, $denovo_line) = @_;
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $standard_line;
    my ($new_alt, $denovo_info_tag, $dq_info_tag, $dngl, $dngt, $dngq)= $self->convert_denovo_line_to_tags($denovo_line, $alt);
    if($denovo_info_tag) {
        $info .= ";$denovo_info_tag";
    }
    if($dq_info_tag) {
        $info .= ";$dq_info_tag";
    }
    $format .= ":DNGL:DNGT:DNGQ";
    for (my $i =0 ; $i < @samples; $i++) {
        $samples[$i] .= ":" . $dngl->[$i] . ":" . $dngt->[$i] . ":" . $dngq->[$i];
    } 
    $alt = $new_alt;
    $standard_line = join("\t", ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples)) . "\n";
    $self->status_message("merged line: $standard_line");
    return $standard_line;
}


sub objectify_denovo_vcf {
    my ($self, $denovo_vcf) = @_;
    my $fh = Genome::Sys->open_file_for_reading($denovo_vcf);
    my %chrs;
    while(my $line = $fh->getline) {
        next if ($line =~ m/^#/);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
        $chrs{$chr}{$pos}=$line;
    }
    return \%chrs;
}

sub convert_denovo_line_to_tags {
    my ($self, $line, $stand_alts) = @_;
    my (@dngt, @dngl, @dngq);
    chomp($line);
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
    my $dq_info_tag;
    for my $info_tag (split(";", $info)) {
        if($info_tag =~m/DQ=/) {
            $dq_info_tag = $info_tag;
        }
    }
    my @alts = split ",", $stand_alts;
    unless($ref eq $stand_alts) {
        unshift @alts, $ref;
    }
    my @denovo_alts = split ",", $alt;
    unless($ref eq $alt) {
        unshift @denovo_alts, $ref;
    }
    my ($gt_idx, $gl_idx, $gq_idx);
    my @format_codes = split ":", $format;
    for (my $i = 0; $i < @format_codes; $i++) {
        if ($format_codes[$i] eq 'GT') {
            $gt_idx=$i;
        }
        if ($format_codes[$i] eq 'GL') {
            $gl_idx = $i;
        }
        if ($format_codes[$i] eq 'GQ') {
            $gq_idx = $i;
        }
    }
    my @denovo_alleles;
    for my $sample (@samples) {
        my @fields = split ":", $sample;
        my $gt = $fields[$gt_idx];
        my ($all1, $all2) = split "/", $gt;
        for my $all ($all1, $all2) {
            unless(grep {/$all/} @alts) {
                push @alts, $all;
            }     
            unless(grep {$_ eq $all} @denovo_alts) { #not in the ALT of the line we were handed
                unless(grep {$_ eq $all} @denovo_alleles) { #not already found by us in another sample or this sample as a homozygote
                        push @denovo_alleles, $all;
                    }
                }     
        }

        my ($all1_idx, $all2_idx);
        for(my $i = 0; $i < @alts; $i++) {
            if($all1 eq $alts[$i]) {
                $all1_idx = $i;
            }
            if($all2 eq $alts[$i]) {
                $all2_idx = $i;
            }
        }
        my $new_gt = "$all1_idx/$all2_idx";
        push @dngt, $new_gt;
        push @dngl, $fields[$gl_idx];
        push @dngq, $fields[$gq_idx];
    }
#    if(@denovo_alleles > 1) {
#        $self->error_message("Two different denovo alleles? $line");
#        die;
#    }
    my $denovo_info_tag = undef;
    if(scalar(@denovo_alleles)>= 1) {
        my @da_indices;
        for my $denovo_allele (@denovo_alleles) {
            my ($index) = grep { $alts[$_] eq $denovo_allele } 0..$#alts;
            push @da_indices, $index;
        }
        $denovo_info_tag =  "DA=" . join(",", @da_indices);
    }
    shift @alts; #throw away ref base, or alt that was same as ref
    $alt = join (",", @alts);
    # called as : my ($new_alt, $denovo_info_tag, \@dngl, \@dngt) = $self->convert_denovo_line_to_tags($denovo_line, $alt);
    return ($alt, $denovo_info_tag, $dq_info_tag, \@dngl, \@dngt, \@dngq);
}



sub fix_denovo_vcf_line {
    my ($self, $line) = @_;
    chomp($line);
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
    my @alts = split ",", $alt;
    unless($ref eq $alt) {
        unshift @alts, $ref;
    }
    my ($gt_idx);
    my @format_codes = split ":", $format;
    for (my $i = 0; $i < @format_codes; $i++) {
        if ($format_codes[$i] eq 'GT') {
            $gt_idx=$i;
        }
    }
    my @denovo_alleles;
    for my $sample (@samples) {
        my @fields = split ":", $sample;
        my $gt = $fields[$gt_idx];
        my ($all1, $all2) = split "/", $gt;
        unless(grep {/$all1/} @alts) {
            push @alts, $all1;
            push @denovo_alleles, $all1;
        }     
        unless(grep {/$all2/} @alts) {
            push @alts, $all2;
            push @denovo_alleles, $all2;
        }     
        my ($all1_idx, $all2_idx);
        for(my $i = 0; $i < @alts; $i++) {
            if($all1 eq $alts[$i]) {
                $all1_idx = $i;
            }
            if($all2 eq $alts[$i]) {
                $all2_idx = $i;
            }
        }
        my $new_gt = "$all1_idx/$all2_idx";
        $fields[$gt_idx]=$new_gt;
        $sample = join(":", @fields);
    }
#    if(@denovo_alleles > 1) {
#        $self->error_message("Two different denovo alleles? $line");
#        die;
#    }
    if(@denovo_alleles>= 1) {
        my @da_indices;
        for my $denovo_allele (@denovo_alleles) {
            my ($index) = grep { $alts[$_] eq $denovo_allele } 0..$#alts;
            push @da_indices, $index;
        }
        my $denovo_info_tag =  "DA=" . join(",", @da_indices);
        $info.= ";$denovo_info_tag";
    }
#    if(@denovo_alleles == 1) {
#        my ($index) = grep { $alts[$_] eq $denovo_alleles[0] } 0..$#alts;
#        my $denovo_info_tag =  "DA=" . $index;
#        $info.= ";$denovo_info_tag";
#    }
    shift @alts; #throw away ref base, or alt that was same as ref
    $alt = join (",", @alts);
    $format =~ s/GT/DNGT/;
    $format = "GT:" . $format;
    $format =~ s/GL/DNGL/;
    $format =~ s/GQ/DNGQ/; 
    for my $sample (@samples) {
        $sample = ".:" . $sample;
    }
    $line = join("\t", ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples));
    return "$line\n";
} 


sub fix_vcf_header {
    my ($self, $standard_vcf) = @_;
    my ($ofh, $temp_output_filename) = Genome::Sys->create_temp_file();

    my $fh = open_vcf_file($standard_vcf);
    while(my $line = $fh->getline) {
        if($line =~ m/^#/) {
            if($line =~m/ID=PS/) {
                $line =~ s/Integer/Float/;
            }   
            if($line =~m/ID=GL/) {
                $line =~s/Unsigned Char, D/Integer,D/;
            }   
            if($line =~m/ID=DS/) {
                $line =~s/Float, D/Float,D/;
            }   
        }
        $ofh->print($line);
    }
    $ofh->close();
    return $temp_output_filename;
}







1;
