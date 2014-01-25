package Genome::Model::Tools::Vcf::VcfToBurdenMatrix;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf "open_vcf_file";
use IO::File;
use FileHandle;
use File::Copy "mv";
use Genome::Sys;
use Carp "croak";

class Genome::Model::Tools::Vcf::VcfToBurdenMatrix {
    is => 'Command::V2',
    has => [
    output_file => {
        is => 'Text',
        is_optional => 0,
        doc => "Output variant matrix format",
        is_input => 1,
        is_output => 1,
    },
    vcf_file => {
        is => 'Text',
        is_optional => 0,
        doc => "Merged Multisample Vcf containing mutations from all samples",
        is_input => 1,
    },
    sample_list_file => {
        is => 'Text',
        is_optional => 1,
        doc => "Limit Samples in the Variant Matrix to Samples Within this File - Sample_Id should be the first column of a tab-delimited file, all other columns are ignored",
    },
    transpose=> {
        is => 'Boolean',
        is_optional => 1,
        doc => "attempt to flip the matrix so that rows  are people, columns are variants, takes more memory",
        default=>0,
        is_input => 1,
    },
    line_buffer_number => {
        is => 'Integer',
        is_optional => 0,
        doc => "Number of lines to buffer before transposing and appending columns to the file. Decrease this number if you are running out of memory.",
        default => 1000,
    },
    vep_annotation_file => {
        is => 'FilesystemPath',
        doc => 'Annotation file to determine which alleles meet requirements',
        is_input => 1,
    },
    burden_test_annotation_file => {
        is => 'Path',
        doc => 'Output annotation file of single or multiple annotation lines per variant',
        is_input => 1,
        is_output => 1,
    },
    consequence_types => {
        is => 'Text',
        doc => 'Which types of alterations to include in the matrix e.g.  ESSENTIAL_SPLICE_SITE,STOP_GAINED,STOP_LOST,NON_SYNONYMOUS_CODING',
        is_many => 1,
        default => ['ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','NON_SYNONYMOUS_CODING'],
    },
    ],
};


sub help_brief {
    "Input Merged Multisample Vcf, Output Variant Matrix for Statistical Genetics Programs"
}


sub help_synopsis {
    <<'HELP';
Input Merged Multisample Vcf, Output Variant Matrix for Statistical Genetics Programs
HELP
}

sub help_detail {
    <<'HELP';
Input Merged Multisample Vcf, Output Variant Matrix for Statistical Genetics Programs
HELP
}

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $vcf_file = $self->vcf_file;
    my $output_file = $self->output_file;
    my $vep_file = $self->vep_annotation_file;
    my $regex_contents = join("|",$self->consequence_types);
    my $type_regex = qr/$regex_contents/;
    #log sample names we're limiting to
    my %samples_to_include;
    if(defined($self->sample_list_file)) {
        my $fh = Genome::Sys->open_file_for_reading($self->sample_list_file);
        while(my $sample_line = $fh->getline ) {
            chomp($sample_line);
            my ($sample_name) = split(/\t/, $sample_line);
            $samples_to_include{$sample_name} = 1;
        }
    }

    #set up anno file
    my $anno_fh = Genome::Sys->open_file_for_reading($vep_file); #this throws on fail so no check
    my ($anno_header, $anno_fetcher) = $self->create_vep_anno_fetcher($anno_fh);

    my $ifh = open_vcf_file($vcf_file);
    my $ofh = Genome::Sys->open_file_for_writing($output_file);
    unless($ofh) {
        die $self->error_message("Unable to open $output_file for writing.");
    }
    if($self->transpose) {
        #want to ensure we have an empty file!
        $ofh->close;
        $ofh = undef;
    }

    my $anno_ofh = IO::File->new($self->burden_test_annotation_file, "w");
    unless($anno_ofh) {
        die $self->error_message("Unable to open " . $self->burden_test_annotation_file. " for writing.");
    }
    $anno_ofh->write(sprintf("#%s\n", join("\t", @$anno_header)));

    my @finished_file;
    my @sample_indexes_to_include;
    while(my $line = $ifh->getline ) {
        chomp($line);
        if ($line =~ m/^\#\#/) {
            next;
        }
        elsif ($line =~ /^\#CHROM/) { #grab sample names off of the header line
            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);
            if(defined($self->sample_list_file)) {
                #initialize which indexes we are keeping out of the sample array
                my $sample_index = 0;
                foreach my $sample (@samples) {
                    if (defined($samples_to_include{$sample})) {
                        push @sample_indexes_to_include, $sample_index;
                    }
                    $sample_index++;
                }
            }
            else {
                @sample_indexes_to_include = 0..scalar(@samples)-1;
            }
            if($self->transpose) {
                push @finished_file, ["VariantId", @samples[@sample_indexes_to_include]];
            }
            else {
                my $header_line = join ("\t", ("VariantId", @samples[@sample_indexes_to_include]));
                $ofh->print("$header_line\n");
            }
            next;
        }
        elsif ($line =~ m/^\#/) { #skip any other commented lines
            next;
        }

        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

        unless($filter =~ m/PASS/i || $filter eq '.') {
            print "Skipping $chr:$pos:$ref/$alt for having filter status of $filter\n";
            next;
        }

        my @sample_list;
        if(defined($self->sample_list_file)) {
            @sample_list = @samples[@sample_indexes_to_include]
        }
        else {
            @sample_list = @samples;
        }

        #parse format line to find out where our pattern of interest is in the ::::: system
        my (@format_fields) = split(/:/, $format);
        my $gt_location; #genotype
        my $ft_location; #filter
        my $count = 0;
        foreach my $format_info (@format_fields) { 
            if ($format_info eq 'GT') {
                $gt_location = $count;
            }
            if ($format_info eq 'FT') {
                $ft_location = $count;
            }
            $count++;
        }
      
        #this file doesn't work if there are unknown genotype locations
        unless ($gt_location == 0) {
            die "Format field doesn't have a GT entry as its first field. No genotypes available or VCF doesn't meet spec for $line\n";
        }
        unless (defined $ft_location) {
            $ft_location = 'NA';
        }
        #grab annotation
        my $annotation_entries = $anno_fetcher->($chr, $pos);
        my ($lines, $anno_line) = $self->format_numeric_output($chr, $pos, $ref, $alt, $id, $gt_location, $ft_location, \@sample_list, $annotation_entries, $type_regex);    
        unless($line && $anno_line) {
            next;
        }
        for my $anno (@$anno_line) {
            print $anno_ofh $anno,"\n";
        }
        if($self->transpose) {
            push @finished_file, @$lines;
            if(scalar @finished_file >= $self->line_buffer_number) {
                #flush buffer
                my $transposed = $self->transpose_row(\@finished_file); #FIXME this requires two copies of the array.
                $self->append_columns_to_file($output_file,$transposed,"\t");
                undef @finished_file;
            }

        } else {
            for my $result_line (@$lines) {
                my $out_line = join("\t", @$result_line);
                $ofh->print("$out_line\n");
            }
        }
    }
    #flush remaining lines
    if($self->transpose) {
        my $transposed = $self->transpose_row(\@finished_file); #FIXME this requires two copies of the array.
        $self->append_columns_to_file($output_file,$transposed,"\t");
    }
    $ofh->close if $ofh;
    return 1;
}

1;

sub find_most_frequent_alleles { 
    my ($self, $chr, $pos, $ref, $alt, $gt_location, $ft_location, $sample_ref) = @_;

    my %alleles_hash;
    my @allele_options;
    foreach my $sample_info (@$sample_ref) {                    
        my (@sample_fields) = split(/:/, $sample_info);
        my $genotype = $sample_fields[$gt_location];
        my $allele1 = my $allele2 = ".";
        ($allele1, $allele2) = split(/\//, $genotype);

        my $filter_status;
        if ($ft_location eq 'NA') {
            $filter_status = 'PASS';
        }
        else {
            $filter_status = $sample_fields[$ft_location];
        }
        $filter_status = defined $filter_status ? $filter_status : 'PASS';

        $sample_info = defined $sample_info ? $sample_info : '.';
        if ($sample_info eq '.') {
        }
        elsif ($filter_status ne 'PASS' && $filter_status ne '.') {
        }
        elsif ($allele1 =~ m/\d+/) {
            $alleles_hash{$allele1}++;
            if ($allele2 =~ m/\d+/) {
                $alleles_hash{$allele2}++;
            }
        }
    }
    @allele_options = (sort { $a <=> $b } keys %alleles_hash);
    my $variant_name;
    my $temp_alt = $alt;
    $temp_alt =~ s/[^A-Za-z0-9]/_/g;    #R itself changes things to a period, but let's try and be specific
    $variant_name = "$chr"."_"."$pos"."_"."$ref"."_"."$temp_alt";

    return ($variant_name, sort { $a <=> $b } keys %alleles_hash);
}


sub format_numeric_output {
    my ($self, $chr, $pos, $ref, $alt, $id, $gt_location, $ft_location, $sample_ref, $annotation_entries, $type_regex) = @_;
    my @samples = @$sample_ref;
    my @alts = split /,/, $alt;
    my @alleles = ($ref, @alts);

    my $anno = $self->get_interesting_variant_per_gene($annotation_entries, $type_regex);
    unless(keys %$anno) {
        return;
    }
    
    #we know we have at least one variant of desired type
    my ($variant_name, @allele_options) = $self->find_most_frequent_alleles($chr, $pos, $ref, $alt, $gt_location, $ft_location, $sample_ref);
    if($id ne '.') {
        $variant_name = $id;
    }

    
    if(@allele_options > 2 && $allele_options[0] != 0) {
        die $self->error_message("Unhandled annotation condition for $variant_name. 3 alleles and the reference is not the major allele.");
    }

    #here, for each gene with a valid variant, create a valid line based on it and an annotation line to go with it.
    my %anno_lines;
    my %most_favored_variant_status;
    for my $index (@allele_options) {
        for my $gene (keys %$anno) {
            #these are in order from most frequent to least frequent
            if($index != 0) {
                my $allele = $alleles[$index];
                if(exists($anno->{$gene}{$allele})) {
                    #so we know that the gene has at least one variant that meets our criteria
                    $most_favored_variant_status{$index}{$gene} = 1;
                    my $cleaned_gene = $gene;
                    $cleaned_gene =~ s/[^A-Za-z0-9]/_/g;
                    $anno_lines{$gene} = "${cleaned_gene}_" . $anno->{$gene}{$allele}->[0]->{_line};
                }
                else {
                    $most_favored_variant_status{$index}{$gene} = 0;
                }
            }
            else {
                $most_favored_variant_status{0}{$gene} = 0;    #set ref to 0. Not true if it's not the most frequent allele
            }
        }
    }
    if($allele_options[0] != 0) {
        $self->warning_message("Reference allele is not the most frequent allele. The annotation for the major allele will be used to determine the effect of the minor (reference) allele. The most frequent allele will be assumed to have no effect on burden.");
        for my $gene (keys %$anno) {
            $most_favored_variant_status{0}{$gene} = $most_favored_variant_status{$allele_options[0]}{$gene}; #swap variant status with the most frequent allele
            $most_favored_variant_status{$allele_options[0]}{$gene} = 0;
        }
    }
    my @return_lines = ();
    my @genes = sort keys %$anno;
    for my $gene (@genes) {
        my $has_variant = 0;
        for my $index (@allele_options) {
            $has_variant ||=  $most_favored_variant_status{$index}{$gene};
        }
        next unless $has_variant;
        my $cleaned_gene = $gene;
        $cleaned_gene =~ s/[^A-Za-z0-9]/_/g;
        my @return_line = ("${cleaned_gene}_$variant_name");
        for my $sample_info (@samples) {
            my (@sample_fields) = split(/:/, $sample_info);
            my $genotype = $sample_fields[$gt_location];
            my $allele1 = my $allele2 = ".";
            ($allele1, $allele2) = split(/\//, $genotype);

            my $filter_status;
            if ($ft_location eq 'NA') {
                $filter_status = 'PASS';
            }
            else {
                $filter_status = $sample_fields[$ft_location];
            }
            $filter_status = defined $filter_status ? $filter_status : 'PASS';

            my $allele_count = 0;
            if ($sample_info eq '.' || $sample_info eq '') {
                $allele_count = '.';
            }
            elsif ($filter_status ne 'PASS' && $filter_status ne '.') {
                $allele_count = '.';
            }
            elsif ($allele1 =~ m/\D+/) {
                $allele_count = '.';
            }
            else {
                for my $allele ($allele1, $allele2) {
                    if(exists($most_favored_variant_status{$allele})) {
                        $allele_count += $most_favored_variant_status{$allele}{$gene};
                    }
                    else {
                        die $self->error_message("Odd allele: $allele for $variant_name with sample $sample_info");
                    }
                }
            }
            push @return_line, $allele_count;
        }
        push @return_lines, \@return_line;
    }
    my @anno_lines = @anno_lines{@genes};
    return (\@return_lines, \@anno_lines);
}


sub transpose_and_print {
    my ($self, $fh, $aoa_ref) = @_;
    my @transposed;
    for my $row (@$aoa_ref) {
        for my $column (0 .. $#{$row}) {
            push(@{$transposed[$column]}, $row->[$column]);
        }
    }

    for my $new_row (@transposed) {
        my $out_line = join("\t", @$new_row);
        $fh->print("$out_line\n");
    }
}

sub transpose_row {
    my ($self, $aoa_ref) = @_;
    my @transposed;
    for my $row (@$aoa_ref) {
        for my $column (0 .. $#{$row}) {
            push(@{$transposed[$column]}, $row->[$column]);
        }
    }
    return \@transposed;
}

sub append_columns_to_file {
    my ($self, $file, $contents, $sep) = @_;

    my $temp_output_file = "$file.append";
    my $ofh = Genome::Sys->open_file_for_writing($temp_output_file);
    #open input file
    if(-s $file) {
        my $ifh = Genome::Sys->open_file_for_reading($file);

        while(my $line = $ifh->getline) {
            chomp $line; #assume we're writing out the newline every single time
            my $append_line = shift @$contents;
            print $ofh join($sep,$line,@$append_line), "\n";
        }
        close($ifh);
    }
    else {
        for my $line (@$contents) {
            print $ofh join($sep, @$line),"\n";
        }
    }

    close($ofh);
    #move the output file to overwrite the input file
    unless(mv($temp_output_file, $file)) {
        die $!;
    }

}

sub create_vep_anno_fetcher {
    my ($self, $fh) = @_;
    unless($fh) {
        croak "Can't create buffered annotation search from undefined file handle\n";
    }
    my $current_entry;
    my @vep_header = qw( Uploaded_variation Location Allele Gene Feature Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation Extra);

    while(my $line = $fh->getline) {
        next if $line =~ /^##/; #skip the meta info
        chomp $line;
        if($line =~ s/^#//) {
            @vep_header = split /\t/, $line;
            last;
        }
        else {
            #first non header line of file. Better go ahead and store it.
            $self->error_message("No header found in VEP file. Using hard-coded header names.");
            my %entry = $self->_make_vep_entry(\@vep_header, $line);
            $current_entry = \%entry;
        }

    }
    my @results;

    my $func = sub {
        my ($chr, $pos) = @_;

        if($current_entry) {
            my ($achr, $apos) = split /:/, $current_entry->{Location};
            if( ($chr eq $achr) && ($pos == $apos) ){
                push @results, $current_entry;
                $current_entry = undef;
            }
        }

        while(my $line = $fh->getline) {
            chomp $line;

            my %entry = $self->_make_vep_entry(\@vep_header,$line);

            my ($achr, $apos) = split /:/, $entry{Location};

            if($achr eq $chr && $apos == $pos) {
                push @results,\%entry;
            }
            elsif(($achr eq $chr && $apos > $pos) || ($achr ne $chr && @results)) {
                #we're done grabbing annotation here because the next entry is further along the genome or at the beginning of the next chromosome
                $current_entry = \%entry;
                my @return_results = @results;
                @results = ();
                return \@return_results;
            }
            else {
                #either apos is less than query or you're not yet on that chromosome
                next;
            }
        }
    };
    return (\@vep_header,$func);
}

sub get_interesting_variant_per_gene {
    my ($self, $annotation_entries, $type_regex) = @_;
    #only have annotation for alts
    my %results_per_gene;
    for my $entry (@$annotation_entries) {
        if($entry->{Consequence} =~ $type_regex) {
            push @{$results_per_gene{$entry->{Gene}}{$entry->{Allele}}}, $entry;  #this shouldn't generate a warning
        }
    }
    return \%results_per_gene;
}

sub _make_vep_entry {
    my ($self, $header, $line) = @_;
    my @fields = split /\t/, $line;
    my %entry;
    @entry{@$header} = @fields;

    #use HGNC name if available
    if(my ($hugo) = $entry{Extra} =~ /HGNC=(\S+)(;|$)/) {
        $entry{Gene} = $hugo;
    }
    $entry{Uploaded_variation} =~ s/[^A-Za-z0-9]/_/g;
    $entry{_line} = join("\t",@entry{@$header});    #reconstruct the line with the new gene name
    return %entry;
}

