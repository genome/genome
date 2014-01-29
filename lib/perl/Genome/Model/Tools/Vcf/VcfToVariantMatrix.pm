package Genome::Model::Tools::Vcf::VcfToVariantMatrix;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf qw(open_vcf_file get_vcf_header);
use IO::File;
use Getopt::Long;
use FileHandle;
use File::Copy "mv";
use Genome::Sys;
use Carp qw/confess/;

class Genome::Model::Tools::Vcf::VcfToVariantMatrix {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            doc => "Output variant matrix format",
            is_input => 1,
            is_output => 1,
            is_optional => 0,
        },
        vcf_file => {
            is => 'Text',
            is_optional => 0,
            doc => "Merged Multisample Vcf containing mutations from all samples",
            is_input => 1,
        },
        positions_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Limit Variant Matrix to Sites - File format chr\\tpos\\tref\\talt",
        },
        bed_roi_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Limit Variant Matrix to Sites Within an ROI - Bed format chr\\tstart\\tstop\\tref\\talt",
        },
        sample_list_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Limit Samples in the Variant Matrix to Samples Within this File - Sample_Id should be the first column of a tab-delimited file, all other columns are ignored",
        },
        project_name => {
            is => 'Text',
            is_optional => 1,
            doc => "Name of the project, will be inserted into output file cell A1",
            default => "Variant_Matrix",
        },
        matrix_genotype_version => {
            is => 'Text',
            is_optional => 1,
            valid_values => [ qw( Bases Numerical) ],
            doc => "Whether or not to output the genotype as the number of non-reference alleles (Numeric) or the actual bases (Bases)",
            default => "Bases",
            is_input => 1,
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
            default => 50000,
        },
        biallelic_only => {
            is => 'Boolean',
            is_optional => 0,
            default => 1,
            doc => "Only include sites that are biallelic in the resulting matrix.",
        },
    ],
    has_transient_optional => [
        _temp_files => {
            is => "ARRAY",
            default => [],
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

###############
sub execute {                               # replace with real execution logic.
    my $self = shift;
    my $vcf_file = $self->vcf_file;
    my $output_file = $self->output_file;
    my $project_name = $self->project_name;
    my $roi_bed = $self->bed_roi_file;
    $DB::single = 1;
    my %sample_name_hash;
    if(defined($self->sample_list_file)) {
        my $sample_list_file = $self->sample_list_file;
        my $sample_list_inFh = Genome::Sys->open_file_for_reading($sample_list_file);
        while(my $sample_line = $sample_list_inFh->getline ) {
            chomp($sample_line);
            my @line_stuff = split(/\t/, $sample_line);
            my $sample_name = $line_stuff[0];
            $sample_name_hash{$sample_name}++;
        }
    }

    my $inFh_positions;
    my %positions_selection_hash;
    if ($self->positions_file) {
        $self->load_positions($self->positions_file, \%positions_selection_hash);
    }

    my ($tfh,$temp_path, $inFh);
    if ($self->bed_roi_file) {
        ## Build temp file for vcf limited to roi ##
        ($tfh,$temp_path) = Genome::Sys->create_temp_file;
        unless($tfh) {
            $self->error_message("Unable to create temporary file $!");
            die;
        }
        $temp_path =~ s/\:/\\\:/g;

        print "Loading Position Restriction File\n";

        my $header = get_vcf_header($vcf_file);
        $tfh->print("$header");
        my $intersect_bed_vcf = `intersectBed -a $vcf_file -b $roi_bed`;
        $tfh->print("$intersect_bed_vcf");
        close ($tfh);
        $inFh = Genome::Sys->open_file_for_reading( $temp_path );
    }
    else {
        $inFh = open_vcf_file($vcf_file);
    }
    print "Loading Genotype Positions from Vcf\n";

    my $fh = IO::File->new($output_file, ">");  # Not using genome sys writing garbage because its a piece of trash
    if($self->transpose) {
        #want to ensure we have an empty file!
        $fh->close;
        $fh = undef;
    }
#########prep done begin main vcf parsing loop############

    my @finished_file;
    my %sample_inclusion_hash;
    while(my $line = $inFh->getline ) {
        chomp($line);
        if ($line =~ m/^\#\#/) {
            next;
        }
        elsif ($line =~ /^\#CHROM/) { #grab sample names off of the header line
            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);
            my @sample_names;
            if(defined($self->sample_list_file)) {
                my $sample_count = 0;
                foreach my $sname (@samples) {
                    if (defined($sample_name_hash{$sname})) {
                        push(@sample_names,$sname);
                        $sample_inclusion_hash{$sample_count}++;
                    }
                    $sample_count++;
                }
            }
            else {
                @sample_names = @samples;
            }
            if($self->transpose) {
                push @finished_file, [$project_name, @sample_names];
            }
            else {
                my $header_line = join ("\t", ($project_name, @sample_names));
                $fh->print("$header_line\n");
            }
            next;
        }
        elsif ($line =~ m/^\#/) { #skip any other commented lines
            next;
        }

        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);
        if ($id || $id ne '.') {
            $id =~ s/,/_/g;
        }

        unless($filter =~ m/PASS/i || $filter eq '.') {
            print "Skipping $chr:$pos:$ref/$alt for having filter status of $filter\n";
            next;
        }

        my @sample_list;
        if(defined($self->sample_list_file)) {
            my $sample_count = 0;
            my @included_samples;
            foreach my $sincluded (@samples) {
                if (defined($sample_inclusion_hash{$sample_count})) {
                    push(@included_samples,$sincluded);
                }
                $sample_count++;
            }
            @sample_list = @included_samples;
        }
        else {
            @sample_list = @samples;
        }

        #parse format line to find out where our pattern of interest is in the ::::: system
        my (@format_fields) = split(/:/, $format);
        my $gt_location; #genotype
        my $dp_location; #depth - this is not filled for vasily ucla files - this is currently unused down below
        my $gq_location; #genotype quality - this is only filled in washu vcf for samtools variants, not varscan and not homo ref - this is currrently unused down below
        my $ft_location; #filter
        my $count = 0;
        foreach my $format_info (@format_fields) {
            if ($format_info eq 'GT') {
                $gt_location = $count;
            }
            elsif ($format_info eq 'DP') {
                $dp_location = $count;
            }
            elsif ($format_info eq 'GQ') {
                $gq_location = $count;
            }
            elsif ($format_info eq 'FT') {
                $ft_location = $count;
            }
            $count++;
        }

        #this file doesn't work if there are unknown genotype locations
        unless ($gt_location || $gt_location == 0) {
            die "Format field doesn't have a GT entry, failed to get genotype for $line\n";
        }
        unless (defined $ft_location) {
            $ft_location = 'NA';
        }
        my $line;
        if($self->matrix_genotype_version=~ m/numerical/i) {
            $line = $self->format_numeric_output($chr, $pos, $id, $ref, $alt, $gt_location, $ft_location, \@sample_list);
        }
        elsif($self->matrix_genotype_version=~ m/bases/i) {
           $line = $self->format_basic_output($chr, $pos, $id, $ref, $alt, $gt_location, $ft_location, \@sample_list);
        }
        else {
            die "Please specify a proper matrix_genotype_version of either \"Bases\" or \"Numerical\"";
        }
        next unless $line;
        if($self->transpose) {
            push @finished_file, $line;
            if(scalar @finished_file >= $self->line_buffer_number) {
                #flush buffer
                my $transposed = $self->transpose_row(\@finished_file); #FIXME this requires two copies of the array.
                $self->append_columns_to_file($transposed,"\t");
                undef @finished_file;
            }

        } else {
            my $out_line = join("\t", @$line);
            $fh->print("$out_line\n");
        }
    }
    #flush remaining lines
    if($self->transpose) {
        my $transposed = $self->transpose_row(\@finished_file); #FIXME this requires two copies of the array.
        $self->append_columns_to_file($transposed,"\t");
        $self->_write_final_output($output_file, "\t");
    }
    $fh->close if $fh;
    return 1;
}

1;

sub load_positions {
    my ($self, $positions_file, $positions_selection_hash) = @_;
    $self->debug_message("Loading Position Restriction File");
    my $inFh_positions = Genome::Sys->open_file_for_reading( $positions_file ) || die "can't open $positions_file\n";
    while(my $line = $inFh_positions->getline ) {
        chomp($line);
        my ($chr, $pos, $ref, $alt) = split(/\t/, $line);
        my $variant_name = "$chr"."_"."$pos"."_"."$ref"."_"."$alt";
        $positions_selection_hash->{$variant_name}++;
    }
}

sub find_unfiltered_alleles {
    my ($self, $chr, $pos, $ref, $alt, $gt_location, $ft_location, $sample_ref) = @_;

    my %alleles_hash;
    for my $sample_info (@$sample_ref) {
        my @sample_fields = split(/:/, $sample_info);

        if ($ft_location ne 'NA' && defined $sample_fields[$ft_location]) {
            next if $sample_fields[$ft_location] ne 'PASS' and $sample_fields[$ft_location] ne '.'
        }
        
        my $genotype = $sample_fields[$gt_location];
        for my $allele (split /\//, $genotype) {
            if($allele =~ m/^\d+$/) {
                $alleles_hash{$allele} = 1;
            }
        }
    }

    my $variant_name = join("_", $chr,$pos,$ref,$alt); 
    $variant_name =~ s/,/_/g;   #remove any commas

    return ($variant_name, sort {$a <=> $b} keys %alleles_hash);
}

sub format_numeric_output {
    my ($self, $chr, $pos, $id, $ref, $alt, $gt_location, $ft_location, $sample_ref) = @_;
    my @samples = @$sample_ref;
    my ($variant_name, @allele_options) = $self->find_unfiltered_alleles($chr, $pos, $ref, $alt, $gt_location, $ft_location, $sample_ref);
    $variant_name .= "_$id" if $id && $id ne '.';
    
    #skip non-biallelic sites
    if($self->biallelic_only && @allele_options != 2) {
        return;
    }

    my $lut = $self->_create_additive_coding_lut(@allele_options);
    
    my @return_line = ($variant_name);
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

        my $allele_count;
        if ($sample_info eq '.' || $sample_info eq '') {
            $allele_count = 'NA';
        }
        elsif ($filter_status ne 'PASS' && $filter_status ne '.') {
            $allele_count = 'NA';
        }
        elsif ($allele1 =~ m/\D+/) {
            $allele_count = 'NA';
        }
        else { 
            $allele_count = $lut->{$allele1}{$allele2};

            unless(defined($allele_count)){
                $allele_count = 'NA';
                print "Couldn't determine allele count for $variant_name with sample info $sample_info\n";
            }
        }
        push @return_line, $allele_count;
    }
    return \@return_line;
}

sub _create_additive_coding_lut {
    my ($self, @alleles) = @_;

    my %coding_lut;
    #enumerate all possible combinations
    #coding number of non-reference alleles
    for my $allele1 (@alleles) {
        for my $allele2 (@alleles) {
            my $total_non_reference = scalar(grep { $_ > 0 } $allele1, $allele2);
            $coding_lut{$allele1}{$allele2} = $total_non_reference;
            $coding_lut{$allele2}{$allele1} = $total_non_reference;
        }
    }
    return \%coding_lut;
}

sub format_basic_output {
    my ($self, $chr, $pos, $id, $ref, $alt, $gt_location, $ft_location, $sample_ref) = @_;
    my @alt_bases = split(/,/, $alt);
    my @allele_option_bases = ($ref, @alt_bases);
    my $variant_name = "$chr"."_"."$pos"."_"."$ref"."_"."$alt";
    $variant_name .= "_$id" if $id && $id ne '.';
    my @return_line = ($variant_name);
    for my $sample_info (@$sample_ref) {
        my (@sample_fields) = split(/:/, $sample_info);
        my $genotype = $sample_fields[$gt_location];
        my ($allele1, $allele2) = split "[/|]", $genotype;

        my $filter_status;
        if ($ft_location eq 'NA') {
            $filter_status = 'PASS';
        }
        else {
            $filter_status = $sample_fields[$ft_location];
        }

        my $allele_type;
        if ($sample_info eq '.' || $sample_info eq '') {
            $allele_type = 'NA';
        }
        elsif ($filter_status ne 'PASS' && $filter_status ne '.') {
            $allele_type = 'NA';
        }
        elsif ($allele1 =~ m/^\D+/) { #if allele1 isn't numerical, then it's not vcf spec and must mean a missing value
            $allele_type = 'NA';
        }
        else { #switch numerical genotypes to the ACTG genotypes
            my $a1 = $allele_option_bases[$allele1];
            my $a2 = $allele_option_bases[$allele2];
            $allele_type = join("",sort($a1,$a2)); #I desire perfection. Also the glm would treat A/G and G/A as two different genotypes.
        }
        push @return_line, $allele_type;
    }
    return \@return_line;
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
    my ($self, $contents, $sep) = @_;

    my ($ofh, $path) = Genome::Sys->create_temp_file;
    push(@{$self->_temp_files}, $path);
    $self->debug_message("Creating temp file #" .scalar(@{$self->_temp_files}). ": $path");
    for my $line (@$contents) {
        $ofh->print(join($sep, @$line)."\n");
    }

    $ofh->close();

}

sub _write_final_output {
    my ($self, $destination, $sep) = @_;

    if (@{$self->_temp_files} == 0) {
        confess "Failed to generate variant matrix: no data stored in temp files!";
    } elsif (@{$self->_temp_files} == 1) {
        #move the output file to overwrite the input file
        unless(mv($self->_temp_files->[0], $destination)) {
            die $!;
        }
    } else {
        my $ofh = Genome::Sys->open_file_for_writing($destination);
        # Hopefully we don't need more than the max # of file descriptors here.
        # I'm comfortable with that assumption for now.
        my @fh = map {Genome::Sys->open_file_for_reading($_)} @{$self->_temp_files};
        print "Merging from " . scalar(@fh) . " temp files.\n";
        while (my @lines = map {$_->getline} @fh) {
            my $undefs = grep {!defined $_} @lines;
            last if ($undefs == @fh);
            if ($undefs) {
                confess "Inconsistent data in temporary files!";
            }
            chomp @lines;
            $ofh->print(join($sep, @lines) . "\n");
        }
    }
}



    #current possible fields, Sept 2011: GT:GQ:DP:BQ:MQ:AD:FA:VAQ:FET
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">
##FORMAT=<ID=BQ,Number=A,Type=Integer,Description="Average Base Quality corresponding to alleles 0/1/2/3... after software and quality filtering">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Average Mapping Quality">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth corresponding to alleles 0/1/2/3... after software and quality filtering">
##FORMAT=<ID=FA,Number=1,Type=Float,Description="Fraction of reads supporting ALT">
##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant Quality">
##FORMAT=<ID=FET,Number=1,Type=String,Description="P-value from Fisher's Exact Test">

