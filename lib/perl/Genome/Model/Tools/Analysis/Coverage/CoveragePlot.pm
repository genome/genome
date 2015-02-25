package Genome::Model::Tools::Analysis::Coverage::CoveragePlot;

# add library path to @INC
use lib (
    'C:\Perl\shared',
    '/gscuser/gchang/usr/bin/perl/site/lib',
    '/gscuser/gchang/usr/bin/perl/shared',
);

use strict;
use Genome;
use IO::File;
use warnings;

use File::Basename;
use File::Spec::Functions; 
use FuncUtils;
use CSlib::Tabular;


class Genome::Model::Tools::Analysis::Coverage::CoveragePlot{
    is => 'Command',
    has => [
        roi_file_path => {
	    is => 'String',
	    is_optional => 0,
            doc => 'Region-of-interest file in the BED format',
        },
	
	output_directory => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'Directory to save output files such as PDF file and tables for plotting, and ref-cov stats files',
	},
	
	bam_files => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'Comma-separated list of bam files to generate coverage plots from, Must specify either --bam_files or --refalign_models',
	},
	
	refalign_models => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'Comma-separated list of refalign model IDs to generate coverage plots from, bam_files will override refalign_models if given',
	},
	
	labels => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'Comma-separated list of the names for each bam (or model), otherwise set such as 1, 2, 3, ... in the order',
	},
	
        min_depth_filters => {
            is => 'String',
            doc => 'Comma-separated list of the minimum depths at each position to consider coverage plot, Must define at least 2 depths',
            default_value => [1, 20, 40],
            is_optional => 1,
        },
	
        wingspan => {
            is => 'Integer',
            doc => 'A base pair wingspan value to add +/- of the input regions, Default is no wingspan',
	    is_optional => 1,
        },
	
	
	# for plotting
        bw => {
            is => 'String',
            doc => 'the smoothing bandwidth to be used in density, Default is "nrd0"',
	    is_optional => 1,
        },
    ]
};

sub help_brief {
    "Sequencing coverage analysis tool to create summary tables and coverage plots by runnning the RefCov tookit"
}

sub help_detail {
    "Sequencing coverage analysis tool to create summary tables and coverage plots by runnning the RefCov tookit"
}



sub execute {
    my $self = shift;
    
    #$DB::deep = 200;		# for debug
    
    
    # defines the folder and file names
    my $dir_refcov = "refcov";
    my $filetab = "table.tsv";
    
    my $script2 = '/gscuser/gchang/usr/tool/sum1stat.pl';
    my $script3 = '/usr/bin/Rscript /gscuser/gchang/usr/tool/stat1plot.R';
    
    
    # gets the arguments
    my $roi_file_path = $self->roi_file_path;
    my $output_directory = $self->output_directory;
    
    my $bam_files = $self->bam_files;
    my $refalign_models = $self->refalign_models;
    
    my $labels = $self->labels;
    my $min_depth_filters = $self->min_depth_filters;
    my $wingspan = $self->wingspan;
    
    my $bw = $self->bw;
    $bw = "nrd0" if ($bw =~ /^\s*$/);
    
    
    # as default
    die "roi-file-path required" unless -e $roi_file_path;
    print STDERR sprintf("%d region(s) of interest from %s\n", `wc -l $roi_file_path | cut -d " " -f 1`, $roi_file_path);
    
    if (-e $output_directory && -d $output_directory)
    {
	print STDERR sprintf("Warning: same directory exists and will overwrite outputs: %s\n", $output_directory);
    }
    
    
    if (defined $min_depth_filters)
    {
	if (ref($min_depth_filters) eq "ARRAY")
	{
	    # with the default
	}
	else
	{
	    $min_depth_filters = [ split(/,/, $min_depth_filters) ];
	}
    }
    else
    {
	#die "min-depth-filters required";
	
	# with the default
	$min_depth_filters = [1, 20, 40];
    }
    
    # adds 1 if not included
    my $found = 0;
    foreach my $d (@$min_depth_filters)
    {
	if ($d == 1)
	{
	    $found = 1;
	}
    }
    
    # checks out for min-depth-filters
    push @$min_depth_filters, 1 unless $found;
    $min_depth_filters = [sort {$a <=> $b} @$min_depth_filters];
    
    die "at least 2 min-depth-filter(s) required" unless @$min_depth_filters > 1;
    print STDERR sprintf("%s as min-depth-filters for coverage plots\n", join(", ", @$min_depth_filters));
    
    if (defined $wingspan)
    {
	die "invalid wingspan: $wingspan" unless $wingspan =~ /^\d+$/;
	print STDERR sprintf("%d-bp wingspan for each region of interest\n", $wingspan);
    }
    
    
    # creates a tabular table
    my $list = CSlib::Tabular->new(
                    # an array of column names used in the tabular data
                    header => [ 'id', 'model', 'bam_path', 'refcov_file' ],
                    
                    # an array of sprintf format of each column in the tabular data
                    format => [ '%s', '%s', '%s', '%s' ],
                );
    
    # prepares the bam list and table
    if (defined $bam_files)
    {
	# if defined, this argument will override the refalign-models parameter
	my @bams = split /\,/, $bam_files;
	
	# adds to the list
	foreach my $b (@bams)
	{
	    die "invalid BAM path: $b" unless -e $b;
	    
	    $list->add([ "", "", $b, "" ]);
	}
    }
    elsif (defined $refalign_models)
    {
	# prepares a bam list from the models
	foreach my $id (split(",", $refalign_models))
	{
            my $model = Genome::Model->get(id =>$id);
	    die "Invalid model ID: $id" unless defined $model;
	    
	    $list->add([ "", $model->id, $model->last_complete_build->merged_alignment_result->bam_path, "" ]);
	    
	    # for debug
	    #print $model->name;
	}
    }
    else
    {
	die "bam-files or refalign-models required";
    }
    
    
    # for the label
    if (defined $labels)
    {
	my @ids = split /,/, $labels;
	
	# if the label number is different with that of BAMs
	die sprintf("%d label(s) required", $list->count_row) unless scalar(@ids) == $list->count_row;
	
	# updates the table
	my %check;
	for (my $i=0; $i<$list->count_row; $i ++)
	{
	    # gets the row hash
	    my $h = $list->row_hash_at($i);
	    
	    # checks if the ID is duplicate
	    if (exists $check{$ids[$i]})
	    {
		die "identical labels found: " . $ids[$i];
	    }
	    else
	    {
		$check{$ids[$i]} ++;
	    }
	    
	    # updates with the ID
	    $h->{id} = $ids[$i];
	    $list->updatehash($i, $h);
	}
    }
    else
    {
	# as default
	# updates the table
	for (my $i=0; $i<$list->count_row; $i ++)
	{
	    # gets the row hash
	    my $h = $list->row_hash_at($i);
	    
	    # updates with the ID
	    $h->{id} = $i;
	    $list->updatehash($i, $h);
	}
    }
    
    
    #create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    if ($tempdir)
    {
	print STDERR "$tempdir temporarily created\n";
    }
    else
    {
        $self->error_message("Unable to create temporary file $!");
	
        die;
    }
    
    # creates directories
    mkdir catfile($tempdir, $dir_refcov);
    
    
    
    # for debug
    #if (1)
    #{
	#$list->data->[0]->[3] = "./refcov-debug/NA12878-CRE-KAPA-refcov-std.tsv";
	#$list->data->[1]->[3] = "./refcov-debug/NA12878-CRE-QXT-refcov-std.tsv";
    #}
    #else
    #{
	# /gscuser/gchang/usr/tool/sum1stat.pl ./refcov-debug 1,40 /tmp/gm-genome_sys-2015-02-05_14_48_50--m60W/anonymous0 nrd0
	# /usr/bin/Rscript /gscuser/gchang/usr/tool/stat1plot.R 1,40 /tmp/gm-genome_sys-2015-02-05_14_48_50--m60W/anonymous0
    #}
    
    # step 1. makes a gmt ref-cov standard run
    my @output1;
    for (my $i=0; $i<$list->count_row; $i ++)
    {
	# gets the row hash
	my $h = $list->row_hash_at($i);
	
	# makes a path to merge readcount results
	my $pathstats = catfile $tempdir, $dir_refcov, sprintf("%s-refcov-std.tsv", $h->{id});
	
	# sets the parameters
	my %params = (
		alignment_file_path => $h->{bam_path},
		alignment_file_format => "bam",
		
		roi_file_path => $roi_file_path,
		roi_file_format => "bed",
		
		min_depth_filter => join(",", @$min_depth_filters),		# no ARRAY
		min_base_quality => 0,
		min_mapping_quality => 0,
		
		stats_file => $pathstats,
		print_headers => 1,
		
		wingspan => $wingspan,
		
		#debug => 1,			# for debug
	    );
	
	# makes a gmt ref-cov standard run
	my $cmd1 = Genome::Model::Tools::RefCov::Standard->create(%params);
	
	if ($cmd1->execute)
	{
	    # updates the table
	    $h->{refcov_file} = catfile $dir_refcov, get_filename($pathstats);
	    $list->updatehash($i, $h);
	}
	else
	{
	    die "gmt ref-cov standard failed for " . $h->{id};
	}
	
	# adds into the array
	push @output1, $pathstats;
    }
    
    # for debug
    #`cp -fr $tempdir/$dir_refcov ./`;
    
    
    # saves the summary table
    print STDERR sprintf("Summary table with %d BAMs saved to %s\n", $list->count_row, $filetab);
    
    my $pathtab = catfile $tempdir, $filetab;
    $list->write($pathtab);
    
    
    # step 2. compiles refcov stats outputs
    my $cmd2 = sprintf "%s %s %s %s %s", $script2, catfile($tempdir, $dir_refcov), join(",", @$min_depth_filters), $tempdir, $bw, join(",", @output1);
    print STDERR sprintf("Compiling RefCov stats files in %s\n%s\n", $dir_refcov, $cmd2);
    system($cmd2);
    
    
    # step 3. creates coverage plots
    my $cmd3 = sprintf "%s %s %s", $script3, join(",", @$min_depth_filters), $tempdir;
    print STDERR sprintf("\nCreating coverage plots\n%s\n", $cmd3);
    system($cmd3);
    
    
    # creats the output directory
    unless (-e $output_directory && -d $output_directory)
    {
        mkdir $output_directory;
	print STDERR sprintf("Output directory created: %s\n", $output_directory);
    }
    
    # copies the output file
    print STDERR "Copying output files\n";
    #Genome::Sys->copy_file("$tempdir/*", './');		# does not work
    `cp -fr $tempdir/* $output_directory/`;
    
    
    # for debug to print out the commands
    #print sprintf("Compiling RefCov stats files in %s\n%s\n", $dir_refcov, $cmd2);
    #print sprintf("\nCreating coverage plots\n%s\n", $cmd3);
    
    print STDERR "Done\n";
    
    return 1;
}

