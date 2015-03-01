package Genome::Model::Tools::Analysis::Coverage::CoveragePlot;

use strict;
use Genome;
use IO::File;
use warnings;
use Carp;

use File::Basename;
use File::Spec::Functions;
use POSIX qw(ceil floor log log10);

use Statistics::R;


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
    "Sequencing coverage analysis tool to create summary tables and coverage plots by runnning the ref-cov tool"
}

sub help_detail {
    "Sequencing coverage analysis tool to create summary tables and coverage plots by runnning the ref-cov tool"
}



sub execute {
    my $self = shift;
    
    #$DB::deep = 200;		# for debug
    
    
    # defines the folder and file names
    my $dir_refcov = "refcov";
    my $filetab = "table.tsv";
    my $script3 = catfile dirname(__FILE__), 'CoveragePlot.R';
    
    
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
    my $list = create_hash(
                    header => [ 'id', 'model', 'bam_path', 'refcov_file' ],
                    format => [ '%s', '%s', '%s', '%s' ],
		    data => [],
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
    # `cp -fr $tempdir/$dir_refcov ./`;
    
    
    # saves the summary table
    print STDERR sprintf("Summary table with %d BAMs saved to %s\n", $list->count_row, $filetab);
    
    my $pathtab = catfile $tempdir, $filetab;
    $list->write($pathtab);
    
    
    # step 2. compiles refcov stats outputs
    #my $cmd2 = sprintf "%s %s %s %s %s", $script2, catfile($tempdir, $dir_refcov), join(",", @$min_depth_filters), $tempdir, $bw, join(",", @output1);
    #print STDERR sprintf("Compiling RefCov stats files in %s\n%s\n", $dir_refcov, $cmd2);
    #system($cmd2);
    my $cmd2 = sprintf "sum1stat: %s %s %s %s %s", catfile($tempdir, $dir_refcov), join(",", @$min_depth_filters), $tempdir, $bw, join(",", @output1);
    print STDERR sprintf("Compiling RefCov stats files in %s\n%s\n", $dir_refcov, $cmd2);
    sum1stat(catfile($tempdir, $dir_refcov), $min_depth_filters, $tempdir, $bw, \@output1, $tempdir);
    
    
    # step 3. creates coverage plots
    my $cmd3 = sprintf "/usr/bin/env Rscript %s %s %s", $script3, join(",", @$min_depth_filters), $tempdir;
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



sub sum1stat {
    my ($dir, $mindepths, $dirout, $bandwidth, $files, $tmp_dir) = @_;
    
    # local variables
    # specifies the columns to retrieve from the refcov statistics files
    my @column = ('name', 'total_ref_bases', 'total_covered_bases', 'ave_cov_depth', 'med_cov_depth');
    my @format = ('%s', '%d', '%d', '%.2f', '%s');
    
    # specifies the column and the min-depth parameter for boxplot statistics
    # (default: 'ave_cov_depth')
    my $statcol = 'ave_cov_depth';
    
    # (default: 1)
    my $statdepth = 1;
    
    
    # as default
    # specifies the directory containing refcov statistics files
    $dir = './refcov' unless defined $dir;
    
    # defines the min-depth filtration for statistics
    $mindepths = [1, 20, 40] unless defined $mindepths;
    
    # defines the output directory
    $dirout = './' unless defined $dirout;
    
    # specifies the smoothing bandwidth to be used in density
    # (default: "nrd0")
    $bandwidth = "nrd0" unless defined $bandwidth;
    
    # specifies file names to compile for statistics
    # (default: empty to read all the files in the directory)
    #@files = ();
    
    # specifies the temporary directory
    $tmp_dir = '/tmp/Statistics-R' . time() unless defined $tmp_dir;
    
    
    # create a communication bridge with R and start R
    # warning message: cannot fetch initial working directory: No such file or directory at /usr/share/perl5/File/Temp.pm line 902
    #my $R = Statistics::R->new(r_bin=>'/usr/bin/R', tmp_dir=>$tmp_dir);
    my $R = Statistics::R->new(r_bin=>'/usr/bin/R');
    $R->startR;
    
    
    # creats a hash for min-depth parameters
    my %filter = map { $_ => 0 } @$mindepths;
    croak "$statdepth min-depth required" unless exists $filter{$statdepth};
    
    # makes the gradient color list
    # orange => [255, 165, 0]
    # yellow => [255, 255, 0]
    # green => [0, 255, 0], green3 => [0, 205, 0]
    # blue => [0, 0, 255],
    # magenta => [255, 0, 255]
    my @colors;
    my $ncolor = scalar(@$mindepths);
    if ($ncolor > 1)
    {
        my $ct = create_hash(color => {});
        
        if ($ncolor == 2)
        {
            @colors = ("orange", "green");
        }
        elsif ($ncolor == 3)
        {
            @colors = ("orange", "yellow", "green");
        }
        elsif ($ncolor < 10)
        {
            my $set = $ct->gradient_3colors([255, 165, 0], [255, 255, 0], [0, 205, 0], $ncolor);
            @colors = map($ct->rgb2hex(@{$set->[$_]}), (0 .. $ncolor - 1));
        }
        elsif ($ncolor < 20)
        {
            my $set;
            my $n = int($ncolor/4) + 1;
            
            while(1)
            {
                $set = $ct->gradient_colors($n, [255, 165, 0], [255, 255, 0], [0, 255, 0], [0, 0, 255]);
                
                if (scalar(@$set) >= $ncolor)
                {
                    last;
                }
                else
                {
                    $n ++;
                }
            }
            
            @colors = map($ct->rgb2hex(@{$set->[$_]}), (0 .. $ncolor - 1));
        }
        else
        {
            my $set;
            my $n = int($ncolor/5) + 1;
            
            while(1)
            {
                $set = $ct->gradient_colors($n, [255, 165, 0], [255, 255, 0], [0, 255, 0], [0, 0, 255], [255, 0, 255]);
                
                if (scalar(@$set) >= $ncolor)
                {
                    last;
                }
                else
                {
                    $n ++;
                }
            }
            
            @colors = map($ct->rgb2hex(@{$set->[$_]}), (0 .. $ncolor - 1));
        }
        
        # for debug
        #printf "\nncolor=%d\n", $ncolor;
        #print_array @colors;
    }
    else
    {
        croak "Too small number for plot: " . $ncolor;
    }
    
    croak "Exception in color number: $ncolor" unless $ncolor == scalar(@colors);
    
    
    # gets all the file names in a directory
    if (@$files > 0)
    {
        foreach my $file (@$files)
        {
            unless (-e $file)
            {
                croak "\nNo file found: $file";
            }
        }
    }
    else
    {
        # gets the file list
        my $dh = DirHandle->new($dir);
        while (defined(my $file = $dh->read))
        {
            if ($file =~ /\.tsv$/i)
            {
                push @$files, catfile($dir, $file);
            }
        }
    }
    
    #printf "\n\n%d files found in %s\n", scalar @$files, $dir;
    #printa @$files;
    
    
    # creates a tabular data object
    my $tabstat = create_hash(
                header => [ 'sample', 'nroi', 'min_depth_filter', 'percent_coverage', 'total_ref_bases', 'total_covered_bases', 'average_depth', 'col' ],      # 'percent_coverage_max', 'percent_coverage_min', 'average_depth_max', 'average_depth_min'
                format => [ '%s', '%d', '%s', '%.2f', '%d', '%d', '%.2f', '%s' ],         # '%.2f', '%.2f', '%.2f', '%.2f'
                data => [],
    );
    
    # creates a tabular data object
    my $tabdist = create_hash(
                    header => [ 'sample', 'min_depth_filter', 'mean', 'median', 'column', 'stats', 'n', 'conf', 'out' ],
                    format => [ '%s', '%s', '%.2f', '%s', '%s', '%s', '%d', '%s', '%s' ],
                    data => [],
                );
    
    # creates a tabular data object
    my $tabden = create_hash(
                    header => [ 'sample', 'i', 'x', 'y', 'bw' ],
                    format => [ '%s', '%d', '%.6f', '%.6e', '%.6f' ],
                    data => [],
                );
    
    
    # compiles the columns
    my %count;
    $count{"nsample"} = 0;      # for the total sample number
    foreach my $path (@$files)
    {
        # gets the file name
        my $file = get_filename($path);
        
        # gets the sample name
        # AML31-KAPA-CRE-stats.tsv
        my ($sample);
        if ($file =~ /^([\w\-]+)\-refcov-std/)
        {
            $sample = "$1";
            
            $count{"nsample"} ++;       # $count{"nfile"}
        }
        else
        {
            croak "Exception in file name: " . $file;
            #$count{sprintf("skip\t%s", $file)} ++;
            
            #next;
        }
        
        
        # reads the refcov statistics table
        my $tab = create_hash(
                    header => undef,
                    format => undef,
                    data => undef,
        );
        $tab->from_file($path, 1);       # with the column names at the first line
        
        printf STDERR "\n  %d ROIs from %s", $tab->count_row, $path;
        
        
        # reads the row in the table
        my %stat;
        my @data = ();
        my $check =  0;
        for (my $i=0; $i<$tab->count_row; $i ++)
        {
            # gets the row hash
            my $h = $tab->row_hash_at($i);    # $tab->row_hash($i) with keys
            
            # filters out given the --min-depth-filter parameter
            my $mindepth;
            if (exists $filter{$h->{min_depth_filter}})
            {
                $mindepth = $h->{min_depth_filter};
            }
            else
            {
                next;
            }
            
            
            # checks an input data
            unless ($check)
            {
                foreach my $c (@column, $statcol)
                {
                    croak "No column found: " . $c unless exists $h->{$c};
                }
                
                # no more checking
                $check = 1;
            }
            
            
            # for the statistics
            if (exists $stat{$sample}->{$mindepth})
            {
                # for the coverage and coverage depth
                $stat{$sample}->{$mindepth}->{total_ref_bases} += $h->{total_ref_bases};
                $stat{$sample}->{$mindepth}->{total_covered_bases} += $h->{total_covered_bases};
                $stat{$sample}->{$mindepth}->{sum_cov_depth} += $h->{ave_cov_depth} * $h->{total_ref_bases};
                
                # if with --print-min-max
                #$stat{$sample}->{$mindepth}->{percent_ref_bases_covered_min} = min2($stat{$sample}->{$mindepth}->{percent_ref_bases_covered_min}, $h->{percent_ref_bases_covered});
                #$stat{$sample}->{$mindepth}->{percent_ref_bases_covered_max} = max2($stat{$sample}->{$mindepth}->{percent_ref_bases_covered_max}, $h->{percent_ref_bases_covered});
                #$stat{$sample}->{$mindepth}->{ave_cov_depth_min} = min2($stat{$sample}->{$mindepth}->{ave_cov_depth_min}, $h->{ave_cov_depth});
                #$stat{$sample}->{$mindepth}->{ave_cov_depth_max} = max2($stat{$sample}->{$mindepth}->{ave_cov_depth_max}, $h->{ave_cov_depth});
            }
            else
            {
                # for the coverage and coverage depth
                $stat{$sample}->{$mindepth}->{total_ref_bases} = $h->{total_ref_bases};
                $stat{$sample}->{$mindepth}->{total_covered_bases} = $h->{total_covered_bases};
                $stat{$sample}->{$mindepth}->{sum_cov_depth} = $h->{ave_cov_depth} * $h->{total_ref_bases};
                
                # if with --print-min-max
                #$stat{$sample}->{$mindepth}->{percent_ref_bases_covered_min} = $h->{percent_ref_bases_covered};
                #$stat{$sample}->{$mindepth}->{percent_ref_bases_covered_max} = $h->{percent_ref_bases_covered};
                #$stat{$sample}->{$mindepth}->{ave_cov_depth_min} = $h->{ave_cov_depth};
                #$stat{$sample}->{$mindepth}->{ave_cov_depth_max} = $h->{ave_cov_depth};
            }
            
            # for the count numbers
            $stat{$sample}->{$mindepth}->{nroi} ++;
            
            # for the boxplot statistics
            if ($mindepth == $statdepth)
            {
                push @data, $h->{$statcol};
            }
            
            # for the count numbers
            $count{$mindepth}->{total_ref_bases} += $h->{total_ref_bases};
            $count{$mindepth}->{total_covered_bases} += $h->{total_covered_bases};
            $count{$mindepth}->{sum_cov_depth} += $h->{ave_cov_depth} * $h->{total_ref_bases};
            
            $count{$mindepth}->{nroi} ++;
            $count{nroi} ++;
        }
        
        
        for (my $i=0; $i<@$mindepths; $i ++)
        {
            my $filter = $mindepths->[$i];
            
            # gets the statistics
            next unless exists $stat{$sample}->{$filter};
            my $s = $stat{$sample}->{$filter};
            
            # calculates the percentage of the reference bases covered
            my $coverage = $s->{total_covered_bases} / $s->{total_ref_bases} * 100;
            
            # calculates the averaged sequencing depth
            my $depth = $s->{sum_cov_depth} / $s->{total_ref_bases};
            
            # calculates the range
            #my $coverage_range = $s->{percent_ref_bases_covered_max} - $s->{percent_ref_bases_covered_min};
            #my $depth_range = $s->{ave_cov_depth_max} - $s->{ave_cov_depth_min};
            
            # adds a data to the end of the tabular data
            $tabstat->add([ $sample, $s->{nroi}, $filter,
                    $coverage, $s->{total_ref_bases}, $s->{total_covered_bases},
                    $depth, $colors[$i] ]);      # $stat{$sample}->{percent_ref_bases_covered_max}, $stat{$sample}->{percent_ref_bases_covered_min}, $stat{$sample}->{ave_cov_depth_max}, $stat{$sample}->{ave_cov_depth_min}
            
            
            if ($filter == $statdepth)
            {
                # calculates the boxplot statistics
                $R->send('d <- c(' . join(",",@data) . ')');
                
                $R->send('m <- median(d)');
                $R->send('b <- boxplot.stats(d)');
                
                # adds to the table
                $tabdist->add([$sample, $statdepth, $depth, get1rvar($R, 'm'), $statcol, join(",", get1rvar($R, 'b$stats')), get1rvar($R, 'b$n'), join(",", get1rvar($R, 'b$conf')), join(",", get1rvar($R, 'b$out'))]);
                
                # computes kernel density estimates
                # known errors:
                # (1) bw.SJ(x, method = "ste") : sample is too sparse to find TD
                # (2) while (fSD(lower) * fSD(upper) > 0) { :   missing value where TRUE/FALSE needed   Calls: density -> density.default -> bw.SJ
                if (is_numeric($bandwidth))
                {
                    $R->send("f <- density(d, bw=$bandwidth, kernel=\"gaussian\", n=512)");
                }
                else
                {
                    $R->send("f <- density(d, bw=\"$bandwidth\", kernel=\"gaussian\", n=512)");
                }
                
                
                # gets the density function
                my @x = get1rvar($R, 'f$x');
                my @y = get1rvar($R, 'f$y');
                my $bw = get1rvar($R, 'f$bw');
                
                # adds to the table
                croak "Exception in density estimates (x, y): $bandwidth" unless scalar(@x) > 0 && scalar(@x) == scalar(@y);
                for (my $i=0; $i<@x; $i ++)
                {
                    $tabden->add([ $sample, $i, $x[$i], $y[$i], $bw ]);
                }
            }
        }
    }
    
    
    # if there is no file
    unless ($count{"nsample"} > 0 && $count{nroi} > 0)
    {
        print STDERR "\n\nNo data found to report: $dir";
        exit
    }
    
    
    # prints out the total statistics
    printf STDERR "\n\n[Summary with %d sample(s)]", scalar @$files;
    printf STDERR "\nnsample\t" . $count{nsample};
    printf STDERR "\nnroi\t" . $count{nroi};
    
    print STDERR "\n\nmin_depth_filter\tpercent_ref_bases_covered\ttotal_ref_bases\ttotal_covered_bases\tave_cov_depth";
    foreach my $filter (sort {$a <=> $b} @$mindepths)
    {
        next unless exists $count{$filter};
        
        my $c = $count{$filter}->{total_covered_bases} / $count{$filter}->{total_ref_bases} * 100;;
        my $d = $count{$filter}->{sum_cov_depth} / $count{$filter}->{total_ref_bases};
        
        printf STDERR "\n%d\t%.2f\t%u\t%u\t%.2f", $filter, $c, $count{$filter}->{total_ref_bases}, $count{$filter}->{total_covered_bases}, $d;
    }
    
    
    # sorts the table
    # NOTE: This will determine the sample order in the plots.
    # as default       $tabdist->sort3("median", C_SORT_NUMD, "mean", C_SORT_NUMD, "sample", C_SORT_STR);
    #$tabdist->{_data} = [ map { $_->[0] }
    #        sort {
    #                    $b->[1] <=> $a->[1]
    #                            ||
    #                    $b->[2] <=> $a->[2]
    #                            ||
    #                    $a->[3] cmp $b->[3]
    #               }
    #        map { [$_, $_->[3], $_->[2], $_->[0]] } @{$tabdist->{_data}} ];
    
    # prints out the boxplot statistics
    #printf "\n\n[Boxplot for %s by %s min-depth: %d samples]\n", $statcol, $statdepth, $tabdist->count_row;
    #$tabdist->print_tab;
    $tabdist->write(catfile($dirout, "boxplot.tab"));
    
    
    # sorts the table
    #$tabstat->sort2("sample", C_SORT_STR, "min_depth_filter", C_SORT_NUM);
    $tabstat->{_data} = [ map { $_->[0] }
                sort {
                            $a->[1] cmp $b->[1]
                                    ||
                            $a->[2] <=> $b->[2]
                       }
                map { [$_, $_->[0], $_->[1]] } @{$tabstat->{_data}} ];
    
    # prints out the tabular table in tab delimited text format
    #printf "\n\n[Sample coverage table: %d samples]\n", $count{"nsample"};
    #$tabstat->print_tab;
    $tabstat->write(catfile($dirout, "coverage.tab"));
    
    
    # sorts the table
    #$tabden->sort2("sample", C_SORT_STR, "i", C_SORT_NUM);
    $tabden->{_data} = [ map { $_->[0] }
                sort {
                            $a->[1] cmp $b->[1]
                                    ||
                            $a->[2] <=> $b->[2]
                       }
                map { [$_, $_->[0], $_->[1]] } @{$tabden->{_data}} ];
    
    # prints out the density function
    #printf "\n\n[Density estimate table]\n";
    #$tabden->print_tab;
    $tabden->write(catfile($dirout, "density.tab"));
    
    
    $R->stopR();
    
    
    # removes the temporary directory
    # `rm -fr $tmp_dir/Statistics-R`;
    
    
    # the end of the main subroutine
}


sub get1rvar {
    # gets the variable as a vector
    my ($R) = shift;
    my ($var) = @_;

    # local variables
    my @arr;
    
    # gets the read
    $R->send('print(' . $var . ')');
    my $print = $R->read;
    
    if (defined $print && $print =~ /^\s*\[1\]\s+/)
    {
        # replaces the multiple numbers
        # [1]    1.0  250.5  500.5  750.5 1000.0
        # [1] 500.5
        # [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54 [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72 [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90 [91]  91  92  93  94  95  96  97  98  99 100  51  50  50  52
        $print =~ s/\s*\[\d+\]\s*/ /g;
        $print =~ s/^\s+//g;
        $print =~ s/\s+$//g;
        
        # gets the value(s)
        @arr = split /\s+/, $print;
    }
    else
    {
        # NULL
        # integer(0)
        # Error in print(x) : object 'x' not found
    }
    
    return (@arr > 1) ? @arr : $arr[0];
}


sub create_hash {
    # creates a hash with the given properties
    my (%property) = @_;
    
    # local variables
    
    # as defaults
    # with no property          croak "property hash required" unless keys(%property) > 0;
    
    
    # with the current package namespace
    my $hash = {};
    bless $hash;
    
    # creates a hash
    foreach my $name (keys %property)
    {
        # assigns the property with its default
        $hash->{"_$name"} = $property{$name};
    }
    
    return $hash;
}


sub get_filename {
    # param: path, return: file name and extension
    # gets the file name
    my $filepath = shift;

    # On Unix returns ("baz", "/foo/bar", ".txt") from "/foo/bar/baz.txt"
    my ($name, $path, $ext) = fileparse($filepath, qr/\.[^.]*/);
    return $name . $ext;
}


sub is_numeric {
    # checks if the value is numeric
    my ($value) = @_;

    if ($value =~ /^\s*$/)
    {
	return 0;
    }
    elsif ($value =~ /\s/)
    {
	return 0;
    }
    else
    {
	if ($value =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/)
	{
	    return 1;
	}
	else
	{
	    return 0;
	}
    }
}


sub count_col {
    # gets the row number
    my ($hsob) = shift;

    # local variables
    
    return (defined $hsob->{_header}) ? scalar(@{$hsob->{_header}}) : 0;
}


sub count_row {
    # gets the column number
    my ($hsob) = shift;

    # local variables
    
    return (defined $hsob->{_data}) ? scalar(@{$hsob->{_data}}) : 0;
}


sub hline {
    # prints out the data line with the given sprint format
    my ($hsob) = shift;
    my ($delimit) = @_;

    # local variables
    $delimit = "\t" unless defined $delimit;
    
    return join($delimit, @{$hsob->{_header}});
}


sub add {
    # adds a data to the end of the tabular data
    my ($hsob) = shift;
    my ($data) = @_;

    # local variables
    
    # the data size
    my $size = scalar(@$data);
    croak sprintf("Mismatched array size of data to add [%d:%d mismatch]\n%s", scalar(@{$hsob->{_header}}), $size, join(" ", @$data)) unless $size == scalar(@{$hsob->{_header}});
    
    # applies the sprintf format to each data element
    # updated on 2012-09-08     my @d = map(sprintf($hsob->{_format}->[$_], $data->[$_]), (0 ... $size - 1));
    my @d;
    for (my $i=0; $i<$size; $i ++)
    {
        unless (defined $data->[$i])
        {
            carp sprintf("Uninitialized value in CSlib::Tabular::add [%d, %s]", $i, $hsob->{_header}->[$i]);
        }
        
        push @d, sprintf($hsob->{_format}->[$i], $data->[$i]);
    }
    
    push @{$hsob->{_data}}, \@d;
}


sub addhash {
    # adds a data to the end of the tabular data with a data hash
    my ($hsob) = shift;
    my ($data, $default) = @_;

    # local variables
    
    # as defaults
    croak "data required" unless defined $data;
    #$default = undef unless defined $default;
    
    # create a data array
    my @d;
    for (my $i=0; $i<@{$hsob->{_header}}; $i ++)
    {
        my $name = $hsob->{_header}->[$i];
        
        if (exists $data->{$name})
        {
            carp "Uninitialized value in $name" unless defined $data->{$name};      # for debug
            
            push @d, sprintf($hsob->{_format}->[$i], $data->{$name});
        }
        else
        {
            push @d, sprintf($hsob->{_format}->[$i], $default);
        }
    }
    
    # adds to the table
    push @{$hsob->{_data}}, \@d;
}


sub updatehash {
    # updates the table with a data hash
    my ($hsob) = shift;
    my ($row, $data, $hash) = @_;

    # local variables
    my $c = 0;
    
    # as defaults
    croak "row required" unless defined $row;
    croak "Out of row range: $row" if ($row < 0 || $hsob->count_row <= $row);
    croak "data required" unless defined $data;
    $hash = $hsob->column_index_hash unless defined $hash;
    
    
    # updates the table
    foreach my $col (keys %$data)
    {
        if (exists $hash->{$col})
        {
            $hsob->{_data}->[$row]->[$hash->{$col}] = sprintf($hsob->{_format}->[$hash->{$col}], $data->{$col});
            
            $c ++;
        }
    }
    
    return $c;
}


sub row_hash_at {
    # gets a row as a row hash for faster speed
    my ($hsob) = shift;
    my ($index) = @_;
    
    # local variables
    
    # as defaults
    # $index for faster speed
    
    # gets a row as an array by the row index number (0-based)
    my $row = $hsob->row($index);
    
    return $hsob->row_hash_by($row, $hsob->{_header});
}


sub write {
    # writes the tabular table to a text file
    my ($hsob) = shift;
    my ($path, $header) = @_;

    # local variables
    
    # as defaults
    $header = 1 unless defined $header;

    # writes into a file
    my $fh = FileHandle->new("> $path");
    croak "Cannot write a file: $path" unless (defined $fh);
    
    print $fh $hsob->text_tab($header);
    
    # prevents the warning message in R
    # saying "incomplete final line found by readTableHeader..."
    print $fh "\n";
    
    $fh->close;
}


sub text {
    # gets a text of the tabular table in a text format
    my ($hsob) = shift;
    my ($delimit, $header) = @_;

    # local variables
    
    # as defaults
    $header = 1 unless defined $header;

    # the header
    my $txt = '';
    $txt = $hsob->hline($delimit) if $header;

    # gets the data line with the given sprint format
    foreach my $d (@{$hsob->{_data}})
    {
        $txt .= "\n" if $header;
        $txt .= $hsob->line($d, $delimit);
        
        $header = 1;
    }
    
    return $txt;
}


sub text_tab {
    # gets a text of the tabular table in tab delimited text format
    my ($hsob) = shift;
    my ($header) = @_;

    # local variables

    return $hsob->text("\t", $header);
}


sub row {
    # gets a row as an array by the row index number (0-based)
    my ($hsob) = shift;
    my ($index) = @_;

    # local variables
    
    croak "Incorrect row index number: $index" if $index < 0 || $hsob->count_row <= $index;
    
    return $hsob->{_data}->[$index];
}


sub row_hash_by {
    # gets a row as a row hash (as a basic method for internal use)
    my ($hsob) = shift;
    my ($data, $keys) = @_;
    
    # local variables
    
    # as defaults
    croak "Data array required" unless defined $data;
    $keys = $hsob->{_header} unless defined $keys;
    croak "Wrong key number" unless $hsob->count_col == @$keys;
    
    
    # creates a row hash
    my %row = map { $keys->[$_] => $data->[$_] } (0 .. @$keys - 1);
    
    return \%row;
}


sub column_index_hash {
    # gets the column index hash
    my ($hsob) = shift;

    # local variables
    my %hash;
    
    for(my $i=0; $i<@{$hsob->{_header}}; $i++)
    {
        croak "The same column name found: " . $hsob->{_header}->[$i] if exists $hash{$hsob->{_header}->[$i]};
        $hash{$hsob->{_header}->[$i]} = $i;
    }
    
    return \%hash;
}


sub line {
    # prints out the data line with the given sprint format
    my ($hsob) = shift;
    my ($data, $delimit) = @_;

    # local variables
    $delimit = "\t" unless defined $delimit;
    
    # the data size
    my $size = scalar(@$data);
    croak "Mismatched array size of data to add" unless $size == scalar(@{$hsob->{_header}});
    
    # applies the sprintf format to each data element
    my @d = map(sprintf($hsob->{_format}->[$_], $data->[$_]), (0 ... $size - 1));
    
    return join($delimit, @d);
}


sub _array_size {
    # returns an array with the given size
    my ($hsob) = shift;
    my ($array, $size) = @_;

    # local variables
    my @new;
    
    foreach my $e (@$array)
    {
        push @new, $e if $size > 0;
        $size --;
    }
    
    return @new;
}


sub from_file {
    # creates a tabular data object from a data file
    my ($hsob) = shift;
    my ($path, $hasheader, $delim, $nrow, $ncol) = @_;

    # local variables
    my $null = "";      # "NA" as the null character
    
    # as a deliminator
    $delim = '\t' unless defined $delim;
    # as a default format
    my $format = '%s';
    # the row and column number
    #$nrow
    #$ncol;
    
    
    # opens the file
    my $fh = FileHandle->new;
    croak "Cannot open a file: $path" unless ($fh->open("< $path"));
    
    # for the column names
    if ($hasheader)
    {
        my $l = $fh->getline;
        chomp $l;
        $l =~ s/\s*$//; # just to make it sure
        $l =~ s/^\s*//;
        
        if ($hasheader == 1)
        {
            # reads the first line
            my @names = split /$delim/, $l;
            @names = $hsob->_array_size(\@names, $ncol) if defined $ncol && $ncol > 0;
            
            if (defined $hsob->{_header})
            {
                # checks the header line
                croak "Different column number" unless scalar(@names) == scalar(@{$hsob->{_header}});
                
                for (my $i=0; $i<@names; $i ++)
                {
                    croak "Different column name: " . $hsob->{_header}->[$i] unless $hsob->{_header}->[$i] eq $names[$i];
                }
            }
            else
            {
                # parses the header line
                for (my $i=0; $i<@names; $i ++)
                {
                    # adds the header and format automatically
                    push @{$hsob->{_header}}, $names[$i];
                    push @{$hsob->{_format}}, $format;
                }
            }
        }
        elsif ($hasheader == 2)
        {
            # skip reading the first line
        }
    }
    
    # add each data line
    my $n = 0;
    my $warn;
    while (my $l = $fh->getline)
    {
        next if $l =~ /^\s*$/;
        
        # for debug     print $l;
        chomp $l;
        #$l =~ s/\s*$//; # just to make it sure
        #$l =~ s/^\s*//;
        
        # parses the data line
        my @fs = split /$delim/, $l;
        @fs = $hsob->_array_size(\@fs, $ncol) if defined $ncol && $ncol > 0;
        
        # checks the data dimension
        # updated on 2012-05-30
        if (scalar @fs > $hsob->count_col)
        {
            # ignores the extra columns
            carp sprintf("Larger column size and ignored (Line %d)\n%s", $n + 1, $l) unless defined $warn;
            $warn = $n;
            
            # resizes the input array
            @fs = $hsob->_array_size(\@fs, $hsob->count_col)
        }
        elsif (scalar @fs < $hsob->count_col)
        {
            # adds empty columns if the column size is smaller
            carp sprintf("Smaller column size and defaulted (Line %d)\n%s", $n + 1, $l) unless defined $warn;
            $warn = $n;
            
            foreach my $i ($#fs + 1 ... $hsob->count_col - 1)
            {
                $fs[$i] = $null;   # with a null character
            }
        }
        
        # adds a data to the end of the tabular data
        # Note: the format of each data element will be reapplied by the given sprintf format
        $hsob->add(\@fs);
        
        $n ++;
        
        if (defined $nrow && $nrow > 0)
        {
            last if $n >= $nrow;
        }
    }
    
    $fh->close;
    
    return $n;
}


sub gradient_2colors {
    # creates a new gradient from the first color to the last color
    my ($hsob) = shift;
    my ($color1, $color2, $ncolor) = @_;
    
    # local variables
    
    # as defaults
    #$color1 = "white" unless defined $color1;
    #$color2 = "black" unless defined $color2;
    $ncolor = 10 unless defined $ncolor;
    
    
    # converts a color name and hexadecimal string to a rgb color
    my ($red1, $green1, $blue1) = @$color1;
    my ($red2, $green2, $blue2) = @$color2;
    
    # makes a gradient color
    my @colors;
    for (my $i=0; $i<$ncolor; $i ++)
    {
        my $red   = $red1 + floor( ( $i / ($ncolor - 1) ) * ( $red2 - $red1 ) );
        my $green = $green1 + floor( ( $i / ($ncolor - 1) ) * ( $green2 - $green1 ) );
        my $blue  = $blue1 + floor( ( $i / ($ncolor - 1) ) * ( $blue2 - $blue1 ) );
        
        push @colors, [$red, $green, $blue];
    }
    
    return \@colors;
}


sub gradient_colors {
    # as a generalized interface with any number of colors
    my ($hsob) = shift;
    my ($nstep, @colors) = @_;
    
    # local variables
    
    # as defaults
    $nstep = 3 unless defined $nstep && $nstep > 2;
    #@colors = ("white", "black") unless @colors > 1;
    
    
    # creates a new gradient between adjacent colors
    my @cols;
    for(my $i=0; $i<scalar(@colors) - 1; $i++)
    {
        my $cs = $hsob->gradient_2colors($colors[$i], $colors[$i+1], $nstep);
        
        if (@cols == 0)
        {
            # adds into the gradient color array
            @cols = @$cs;
        }
        else
        {
            # joins the two color arrays
            shift @$cs;
            push @cols, @$cs;
        }
    }
    
    #printl scalar(@cols);      # for debug
    return \@cols;
}


sub gradient_3colors {
    # creates a new gradient with 3 colors
    my ($hsob) = shift;
    my ($color1, $color2, $color3, $ncolor) = @_;
    
    # local variables
    
    # calculates the color numbers
    my ($n1, $n2) = ($ncolor % 2 == 0) ? ($ncolor / 2, $ncolor / 2 + 1) : (ceil($ncolor / 2), ceil($ncolor / 2));
    
    # creates a new gradient from the first color to the last color
    my $col1 = $hsob->gradient_2colors($color1, $color2, $n1);
    my $col2 = $hsob->gradient_2colors($color2, $color3, $n2);
    
    # joins the two color arrays
    shift @$col2;
    push @$col1, @$col2;
    
    #printl scalar(@$col1);      # for debug
    return $col1;
}


sub rgb2hex {
    # converts a rgb color to a hexadecimal string
    my ($hsob) = shift;
    my ($red, $green, $blue) = @_;
    
    # local variables
    
    # as defaults
    
    
    # converts a decimal color code into a hex color code
    my $hexcolor = sprintf ("#%02lx%02lx%02lx", $red, $green, $blue);
    
    return uc $hexcolor;
}

