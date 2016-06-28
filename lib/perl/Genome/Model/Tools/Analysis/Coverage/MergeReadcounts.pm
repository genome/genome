package Genome::Model::Tools::Analysis::Coverage::MergeReadcounts;
use strict;
use Genome;
use IO::File;
use warnings;

use List::MoreUtils qw(each_array all pairwise);


class Genome::Model::Tools::Analysis::Coverage::MergeReadcounts{
    is => 'Command::V2',
    has => [
	bam_files => {
	    is => 'String',
	    is_optional => 0,
	    is_many => 1,
	    doc => 'comma-separated list of bam files to grab readcounts from, Output columns will be appended in this order',
	},
	
	variant_files => {
	    is => 'String',
	    is_optional => 0,
	    is_many => 1,
	    doc => 'comma-separated list of text files containing snvs in annotation format (1-based, first 5-cols = [chr, st, sp, ref, var]), all variant sites will be combined to a list',
	},
	
	variant_sources => {
	    is => 'String',
	    is_optional => 1,
	    is_many => 1,
	    doc => 'comma-separated list of the name for each variant_file, used in the source column (1-based index number as default such as 1, 2, 3, ...)',
	},
	
	output_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output file will be indentical to the input file with readcounts appended as the last two columns',
	},
	
	genome_build => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'takes either a string describing the genome build (one of 36, 37, mm9, mus37) or a path to the genome fasta file',
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
	    is => 'Integer',
	    is_optional => 1,
	    doc => 'minimum depth required for a site to be reported',
	},
	
	max_depth => {
	    is => 'Integer',
	    is_optional => 1,
	    doc => 'maximum depth allowed for a site to be reported',
	},
	
	min_vaf => {
	    is => 'Integer',
	    is_optional => 1,
	    doc => 'minimum variant allele frequency required for a site to be reported (0-100)',
	},
	
	max_vaf => { 
	    is => 'Integer',
	    is_optional => 1,
	    doc => 'maximum variant allele frequency allowed for a site to be reported (0-100)',
	},
	
	indel_size_limit => {
	    is => 'Integer',
	    is_optional => 1,
	    doc => 'maximum indel size to grab readcounts for. (The larger the indel, the more skewed the readcounts due to mapping problems)',
	    default => 4,
	},
	
	header_prefixes => { 
	    is => 'String',
	    is_optional => 1,
	    is_many => 1,
	    doc => 'Comma-separated list - if the file has a header, three column titles get added for each bam ("ref_count","var_count","VAF"). This specifies a prefix for those columns. (i.e.  "Normal" will lead to "Normal_ref_count","Normal_var_count","Normal_VAF"). As a default, the ascending integers in the same BAM order, e.g. 1, 2, 3, ...',
	},
	
	per_library => {
	    is => 'Boolean',
	    is_optional => 1,
	    doc => 'whether or not to report counts on a per-library basis',
	    default => 0,
	},
	
	bam_readcount_version => {
	    is => 'String',
	    doc => 'version of bam-readcount to use',
	    is_optional => 1,
	},
	
    ]
};

sub help_brief {
    "get readcounts. make pretty. create anno file by merging anno files"
}

sub help_detail {
    "get readcounts. make pretty. create anno file by merging anno files"
}



sub execute {
    my $self = shift;
    my @bam_files = $self->bam_files;
    my @variant_files = $self->variant_files;
    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $min_quality_score = $self->min_quality_score;

    my $min_vaf = $self->min_vaf;
    my $max_vaf = $self->max_vaf;
    my $min_depth = $self->min_depth;
    my $max_depth = $self->max_depth;
    my $indel_size_limit = $self->indel_size_limit;

    my $chrom = $self->chrom;
    my @header_prefixes;
    my @sources;
    
    if (defined($self->variant_sources))
    {
	@sources = $self->variant_sources;
	
	die "Source and variant file do not match" unless scalar(@variant_files) == scalar(@sources);
    }
    else
    {
	@sources = (1 ... scalar(@variant_files));
    }
    
    if (defined($self->header_prefixes))
    {
        @header_prefixes = $self->header_prefixes;
    }

    # checks the BAM and header prefixes in 1:1 relationship
    if (@header_prefixes > 0)
    {
	unless (scalar(@header_prefixes) == scalar(@bam_files))
	{
	    die "Header and variant file do not match";
	}
    }
    else
    {
	# as a default
	@header_prefixes = (1 ... scalar(@header_prefixes));
    }
    
    my $fasta;
    if ($genome_build eq "36") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif( ($genome_build eq "37") || ($genome_build eq "37lite") ){
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif ($genome_build eq "mus37") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-mouse-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } elsif ($genome_build eq "mm9") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "UCSC-mouse-buildmm9");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } elsif (-e $genome_build ) {
        $fasta = $genome_build;
    } else {
        die ("invalid genome build or fasta path: $genome_build\n");
    }
    

    # The following defaults are for the AddReadcounts module.
    #set some defaults for undefined parameters
    #unless(defined($min_depth)){
    #    $min_depth = "NA";
    #}
    #unless(defined($max_depth)){
    #    $max_depth = "NA";
    #}
    #unless(defined($min_vaf)){
    #    $min_vaf = "NA";
    #}
    #unless(defined($max_vaf)){
    #    $max_vaf = "NA";
    #}
    #unless(defined($chrom)){
    #    $chrom = "all";
    #}



    #create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    unless($tempdir) {
        die("Unable to create temporary file $!");
    }
    
    
    # reads the variant files and combine
    # for a hash with the keys chr, start, stop, ref, var, type
    my (%hash, @header);

    my $ea = each_array(@variant_files, @sources);
    while (my ($variant, $source) = $ea->())
    {
	# creates a file handler
	my $fh = FileHandle->new;
	die "Cannot open a file: " . $variant unless $fh->open("< $variant");
	
	# for the header line
	my $line = $fh->getline;
	chomp $line;
	if (scalar(@header) == 0)
	{
	    @header = split /\t/, $line;
	}
	else
	{
	    my @new = split /\t/, $line;
	    
	    die "Inconsistent table format found in the header line: " . $variant unless scalar(@header) == scalar(@new);
	    
	    # double-checks whether they share the same column header
	    unless (all { $_ } pairwise { $a eq $b } @header, @new)
	    {
		die "Inconsistent table format found in the header line: " . $variant;
	    }
	}
	
	while (my $i = $fh->getline)
	{
	    next if $i =~ /^\s*$/;
	    chomp $i;
	    
	    # gets the field values
	    my ($chr, $start, $stop, $ref, $var, $type, @vs) = split /\t/, $i;
	    
	    # The following line was to chop off the read count columns assuming about the 25 extra columns,
	    # but inactivated in order to keep the original header line.
	    # @vs = splice @vs, 0, (25 - 6);   # prints only 25 fields
	    
	    # creates a hash with the field values
	    my %h = ("chr" => $chr, "start" => $start, "stop" => $stop, "ref" => $ref, "var" => $var, "type" => $type, "extra" => \@vs,
			"source" => {$source => 1} );
	    
	    # creates a mutation site ID
	    my $id = sprintf "%s:%u:%u:%s:%s:%s", $chr, $start, $stop, $ref, $var, $type;
	    
	    # in theory	die "Identical mutation site found: $id" if exists $hash{$id};
	    if (exists $hash{$id})
	    {
		$hash{$id}->{source}->{$source} ++;
	    }
	    else
	    {
		$hash{$id} = \%h;
	    }
	}
	
	$fh->close;
    }
    
    # sorts the mutation sites by their position
    my @sort = map { $_->[0] }
           sort { $a->[1] cmp $b->[1]
                           ||
                    $a->[2] <=> $b->[2]
                           ||
                    $a->[3] <=> $b->[3]
                           ||
                    $b->[4] cmp $a->[4]
           } map { [$_, $hash{$_}->{chr}, $hash{$_}->{start}, $hash{$_}->{stop}, $hash{$_}->{type}] } keys %hash;
    
    # prints out the progress
    $self->status_message("%d mutation sites combined from %d annotation files\n", scalar(@sort), scalar(@variant_files));
    
    
    # writes the combined mutation sites in the annotation format
    my $pathtab = "$tempdir/snvs.indels.annotated-combine";
    my $fh = FileHandle->new("> $pathtab");
    die "Cannot write a file: " . $pathtab unless defined $fh;
    
    # The following line was to chop off the read count columns assuming about the 25 extra columns,
    # but inactivated in order to keep the original header line.
    # @header = splice @header, 0, 25;   # prints only 25 columns in the header line
    
    # writes the header line
    print $fh join("\t", @header, "source") . "\n";
    
    # writes the mutation sites in the sorted order
    foreach my $id (@sort)
    {
	# gets the source value
	my $source = join(",", sort keys %{$hash{$id}->{source}});
	
	print $fh sprintf("%s\t%u\t%u\t%s\t%s\t%s\t%s\n", $hash{$id}->{chr}, $hash{$id}->{start}, $hash{$id}->{stop}, $hash{$id}->{ref}, $hash{$id}->{var}, $hash{$id}->{type}, join("\t", @{$hash{$id}->{extra}}, $source));
    }
    
    $fh->close;
    
    
    # makes a path to merge readcount results
    my $pathmerge = "$tempdir/snvs.indels.annotated-merge";
    
    # makes an AddReadcounts run
    my %params = (
        bam_files => join(",", @bam_files),
        output_file =>  $pathmerge,
        variant_file => $pathtab,
        genome_build => $genome_build, 
        header_prefixes => join(",", @header_prefixes),
        chrom => $chrom,
        min_depth  => $min_depth,
        max_depth => $max_depth,
        min_vaf => $min_vaf,
        max_vaf => $max_vaf,
        indel_size_limit => $indel_size_limit,
        min_quality_score => $min_quality_score,
        per_library => $self->per_library,
        bam_readcount_version => $self->bam_readcount_version,
        );

    if($self->bam_readcount_version){
        $params{bam_readcount_version} = $self->bam_readcount_version;
    }
    my $cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(%params);
    unless ($cmd->execute) {
	die "add-readcounts failed";
    }
    
    # copies the output file
    # copy($pathmerge, $output_file) or die "Copy failed: $!";
    `cp $pathmerge $output_file`;
    
    return 1;
}
