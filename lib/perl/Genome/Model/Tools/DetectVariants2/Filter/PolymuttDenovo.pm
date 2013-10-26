package Genome::Model::Tools::DetectVariants2::Filter::PolymuttDenovo;

use strict;
use warnings;
use Data::Dumper;
use Genome;
use Genome::Utility::Vcf "open_vcf_file";
use Genome::Info::IUB;
use POSIX;
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::DetectVariants2::Filter::PolymuttDenovo {
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has_optional_input => [
    min_read_qual=> {
        doc=>'the lowest quality reads to use in the calculation. default q20',
        default=>"20",
    },
    min_unaffected_pvalue=> {
        doc=>"the minimum binomial test result from unaffected members to pass through the filter",
        default=>"1.0e-4",
    },
    ],
    doc => "A binomial filter for polymutt denovo output",
};

sub help_detail {
    "A binomial filter for polymutt denovo output"
}

sub _variant_type { 'snvs' };

sub _filter_name { 'PolymuttDenovo' };

sub _filter_variants {
    my $self=shift;
    my $vcf = $self->input_directory . "/snvs.vcf.gz"; # TODO this should probably just operate on snvs.vcf.gz and only filter denovo sites (info field?)
    my $output_file = $self->_temp_staging_directory. "/snvs.vcf.gz";

    my @alignment_results = $self->alignment_results;
    $self->error_message("Alignment Results found: " . scalar(@alignment_results));

    my $sites_file = Genome::Sys->create_temp_file_path();
    my $cat_cmd = "cat";
    my $vcf_fh = open_vcf_file($vcf);
    if(Genome::Sys->file_is_gzipped($vcf)) {
        $cat_cmd = "zcat";
    }

    my $sites_cmd = qq/ $cat_cmd $vcf | grep -v "^#" | grep "DA" | awk '{OFS="\t"}{print \$1, \$2, \$2}' > $sites_file/;

    unless(-s $sites_file) {
        print STDERR "running $sites_cmd\n";
        `$sites_cmd`;
    }
    unless(-s $sites_file) {
        $self->status_message("No denovo sites found to filter. this filter is a no-op, copying file over.");
        `cp $vcf $output_file`;
        return 1;
    }


    # Sort the alignment results in a sane way
    my $header_cmd = qq/ $cat_cmd $vcf | grep "^#CHR"/;
    my $header_line= `$header_cmd`;
    unless ($header_line) {
        die $self->error_message("Could not get the header line with subject names using $header_cmd");
    }
    @alignment_results = $self->sort_alignment_results_by_header($header_line, @alignment_results);

    ###prepare readcounts
    my $ref_fasta=$self->reference_sequence_input;
    my $readcount_file_ref = $self->prepare_readcount_files($ref_fasta, $sites_file, \@alignment_results);

    ###prepare R input
    my $r_input_path = $self->prepare_r_input($vcf_fh, $readcount_file_ref);
    ### RUN R
    my ($r_fh, $r_script_path) = Genome::Sys->create_temp_file();
    my ($r_output) = Genome::Sys->create_temp_file_path();
    $r_fh->print($self->r_code($r_input_path, $r_output));
    $r_fh->close;
    $self->status_message(join('', `R --vanilla < $r_script_path 2>&1`));
    ###output new VCF
    my %ped_hash = $self->make_trios($self->pedigree_file_path, $header_line);
    $self->output_passing_vcf($vcf, $r_output, $self->min_unaffected_pvalue, $output_file, \%ped_hash);
    return 1;
}
sub make_trios {
    my ($self, $ped_file, $header_line) = @_;
    chomp($header_line);
    my @fields = split "\t", $header_line;
    splice(@fields, 0, 9);
    my %person_to_index;
    my %ped_hash;
    for (my $i =0; $i < scalar(@fields); $i++) {
        $person_to_index{$fields[$i]}=$i;
    }
    my $fh = IO::File->new($ped_file);
    while(my $line = $fh->getline) {
        my ($family, $individual, $dad, $mom, $sex, $glf, $affected) = split "\t", $line;
        if($dad && $mom) {
            my $child_idx = $person_to_index{$individual};
            $ped_hash{$child_idx}{'dad_index'}=$person_to_index{$dad};
            $ped_hash{$child_idx}{'mom_index'}=$person_to_index{$mom};
        }
        elsif($dad || $mom) {
            $self->error_message("Previously we discarded incomplete peds upstream");
            $self->error_message("If this has changed, this filter needs to change.");
            die();
        }
    }
    return %ped_hash;
}






sub find_indexes_of_unaffecteds {
    my ($self, $header_line, $ped_file) = @_;
    my $ped_fh = IO::File->new($ped_file);
    my @unaffected_indices;
    my %status;
    while(my $line = $ped_fh->getline) {
        chomp($line);
        my @fields = split "\t", $line;
        my $sample_name = $fields[1];
        my $affectation_status = $fields[-1];
        if(($affectation_status ne 'U') && ($affectation_status ne 'A')) {
            $self->error_message("Ped line: $line does not end in U or A. Unable to parse ped to figure out who is unaffected.");
            die;
        }
        $status{$sample_name}=$affectation_status;
    }
    $ped_fh->close();
    my @header = split"\t", $header_line;
    splice(@header, 0, 9);
    for(my $i = 0; $i < scalar(@header); $i++) {
        if(exists($status{$header[$i]})) {
            if($status{$header[$i]} eq 'U') {
                push @unaffected_indices, $i;
            }
        }
    }
    return \@unaffected_indices;
}

sub output_passing_vcf {
    my($self, $input_filename, $r_output, $pvalue_cutoff, $output_filename, $ped_hashref) = @_;
    my $pvalues = Genome::Sys->open_file_for_reading($r_output);
    my %pass_hash;    

    while(my $pvalues_line = $pvalues->getline) { #parse pvalue returns and figure out who is gonna pass
        my @fields = split"\t", $pvalues_line;
        my $pvalues_chr = $fields[0];
        my $pvalues_pos = $fields[1];
        splice(@fields,0,2);
        my $pass = "PASS";
        for my $individual (keys %{$ped_hashref}) {
            my $dad_idx = $ped_hashref->{$individual}->{'dad_index'};
            my $mom_idx = $ped_hashref->{$individual}->{'mom_index'};
            for my $parent ($dad_idx, $mom_idx) {
                my $unaff_pval = $fields[$parent];
                if($unaff_pval < $pvalue_cutoff) {
                    $pass=$self->_filter_name;
                }
            }
            $pass_hash{$pvalues_chr}{$pvalues_pos}{$individual}=$pass;
        }
    }
    my $input_file = Genome::Sys->open_gzip_file_for_reading($input_filename);
    my $output_file = Genome::Sys->open_gzip_file_for_writing($output_filename);
    my $header_printed=0;
    my @filtered_already=`zcat $input_filename | grep "DNFT"`; #we got hacks
    if(@filtered_already) {
        $header_printed=1;
    }    
    while(my $line = $input_file->getline) { #stream through whole file and only touch the members of the hash we filled out above 

        if($line =~m/^#/) {
            if($line=~m/^##FORMAT/ && !$header_printed) {
                $output_file->print(qq|##FORMAT=<ID=DNFT,Number=1,Type=String,Description="Denovo Filter Status">\n|);
                $header_printed=1;
            }
            $output_file->print($line);
            next;  #don't try to filter if this was a header line, asshole
        }
        chomp($line);
        my ($vcf_chr, $vcf_pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split("\t", $line);
        if(exists($pass_hash{$vcf_chr}{$vcf_pos})){
            my $DNFT_idx = $self->find_format_field($format, "DNFT");
            unless(defined($DNFT_idx)) {
                $self->status_message("appending DNFT");
                $format .= ":DNFT";
            }
            for (my $i=0; $i < scalar(@samples); $i++) {
                my $filter_status = $pass_hash{$vcf_chr}{$vcf_pos}{$i} || ".";
                if(defined($DNFT_idx)) {
                    my @sample_fields = split ":", $samples[$i];
                    $sample_fields[$DNFT_idx]=$filter_status;
                    $samples[$i]=join(":", @sample_fields);
                }
                else {
                    $samples[$i] .= ":$filter_status";
                }
            }
            $line = join "\t", ($vcf_chr, $vcf_pos, $id, $ref, $alt, $qual, $filter, $info,$format, @samples);
        }
        $output_file->print("$line\n");

    }
    return 1;
}

sub prepare_readcount_files {
    my $self = shift;
    my $ref_fasta = shift;
    my $sites_file = shift;
    my $alignment_results_ref = shift;
    my @readcount_files;
    my $qual = $self->min_read_qual;
    for my $alignment_result(@$alignment_results_ref) {
        my $sample_name = $self->find_sample_name_for_alignment_result($alignment_result);
        my $readcount_out = Genome::Sys->create_temp_file_path($sample_name . ".readcount.output");
        my $bam = $alignment_result->merged_alignment_bam_path;
        unless (-f $bam) {
            die "merged_alignment_bam_path does not exist: $bam";
        }
        push @readcount_files, $readcount_out;

        #my $readcount_cmd = "bam-readcount -q $qual -f $ref_fasta -l $sites_file $bam > $readcount_out";
        unless(-s $readcount_out) {
            print STDERR "running bam-readcount";
            my $rv = Genome::Model::Tools::Sam::Readcount->execute(
                bam_file => $bam,
                minimum_mapping_quality => $qual,
                output_file => $readcount_out,
                reference_fasta => $ref_fasta,
                region_list => $sites_file,
            );
            unless($rv) {
                $self->error_message("bam-readcount failed");
                return;
            }
        }
    }
    return \@readcount_files;
}

sub prepare_r_input {
    my ($self, $vcf_fh, $readcount_files, $ped_hashref) = @_;
    my ($r_input_fh,$r_input_path) = Genome::Sys->create_temp_file();

    while (my $line = $vcf_fh->getline) {
        next if ($line =~ m/^#/);
        next unless($line =~m/DNGT/); #should be DA?  then we could remove "if ($novel)" below I imagine
        chomp($line);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
        my @alts = split ",", $alt;
        my @possible_alleles = ($ref, @alts);
        my @novel = ();
        my @info_fields = split";", $info; 
        for my $info_field (@info_fields) {
            if($info_field =~m/^DA=/) {
                my ($novel_string) = ($info_field =~m/^DA=(\d+)/);
                my @novel_idx = split ",", $novel_string;
                for my $novel_i (@novel_idx) {
                    push @novel, $possible_alleles[$novel_i];
                }
                $self->status_message("Selecting bases @novel as novel allele for line: $line\n");
            }
        }
        my $DNGT_idx=$self->find_format_field($format, "DNGT");
        unless(defined($DNGT_idx)) {
            $self->error_message("Need DNGT field to proceed. couldn't find one in $line");
            die;
        }
        if(scalar(@novel) >= 1) {
            $r_input_fh->print("$chr\t$pos");
            for (my $i=0; $i < scalar(@samples); $i++) {    
                my $readcount_file = $readcount_files->[$i];
                chomp(my $readcount_line = `grep "^$chr\t$pos" $readcount_file`);
                my ($chr, $pos, $rc_ref, $depth, @fields) = split "\t", $readcount_line;
                my $prob;
                my (@format_fields) = split ":", $samples[$i];
                my $dngt = $format_fields[$DNGT_idx];
                my @ref_indices = split "[|/]", $dngt;
                my @ref_bases = map { $possible_alleles[$_] } @ref_indices;
                if($self->contains(\@novel, \@ref_bases)) {
                    $prob = .5;
                }
                else {
                    $prob = .01;
                }
                my $ref_depth=0;
                my $var_depth=0;
                if($readcount_line) {
                    for my $fields (@fields) {
                        my ($base, $depth, @rest)  = split ":", $fields;
                        next if($base =~m/\+/);
                        if(grep /$base/,@novel) {
                            $var_depth+=$depth;
                        }
                        else {
                            $ref_depth+=$depth;
                        }
                    }
                }
                $r_input_fh->print("\t$var_depth\t$ref_depth\t$prob");
            } #end adding all people for one variant onto one line (3 values per person)
            $r_input_fh->print("\n");
        }  #move to new variant line
    }
    $r_input_fh->close;
    return $r_input_path;
}

sub contains {
    my ($self, $array1, $array2) = @_;
    my %hash;
    for my $value (@$array1) {
        $hash{$value}=1;
    }
    for my $value (@$array2) {
        if(exists($hash{$value})) {
            return 1;
        }
    }
    return 0;
}



sub r_code {
    my $self = shift;
    my $input_name = shift;
    my $output_name = shift;
    return <<EOS
options(stringsAsFactors=FALSE);
mytable=read.table("$input_name");
num_indices=(ncol(mytable)-2)/3
pvalue_matrix=matrix(nrow=nrow(mytable), ncol=num_indices);
indices = seq(from=3,to=num_indices*3, by=3)
for (i in 1:nrow(mytable)) {
    for (j in indices) {
        pvalue_matrix[i,(j/3)]=binom.test(as.vector(as.matrix(mytable[i,j:(j+1)])), 0, mytable[i,(j+2)], "t", .95)\$p.value;
    }
}
chr_pos_matrix=cbind(mytable\$V1,mytable\$V2, pvalue_matrix);
write(t(chr_pos_matrix), file="$output_name", ncolumns=(num_indices+2), sep="\t");
EOS
}

sub find_format_field {
    my ($self, $format, $field) = @_;
    my $idx;
    my @format_fields = split ":", $format;
    for (my $i = 0; $i < scalar(@format_fields); $i++) {
        my $format_tag = $format_fields[$i];
        if($format_tag eq $field) {
            $idx = $i;
        }
    }
    unless(defined($idx)) {
        $self->error_message("unable to find $field in format");
        return undef;
    }
    return $idx;
}



sub _generate_standard_files {
    return 1;
}

1;
