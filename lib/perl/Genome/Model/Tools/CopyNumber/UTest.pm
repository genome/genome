package Genome::Model::Tools::CopyNumber::UTest;

# this is different from backup050610, since we do gene hit check first.

use strict;
use Genome;
use Cwd 'abs_path';
use IO::File;
use Getopt::Long;
use Statistics::R;
use File::Temp;
use DBI;
use Genome::Utility::HugoGene::HugoGeneMethods;
require Genome::Sys;

class Genome::Model::Tools::CopyNumber::UTest {
    is => 'Command',
    has => [
    output_file => {
        is => 'String',
        is_optional => 1,
        doc => 'File name contain the U test P value in the form of "name,chromosome,start-position,end-position,sliding-window,tumor-p-value,normal-p-value,tumor-normal-p-value,tumor_mean,normal_mean,gene_hit,plot-graph. Should include the full directory."',
    },
    input_file => {
    	is => 'String',
    	is_optional => 1,
    	doc => 'File name in the form of "name,tumor-bam-file,normal-bam-file,chromosome,start-position,end-position". Leave it blank if any item does not apply but please do not skip comma. No header needed. Should include the full directory.',
    },
    graph_directory => {
    	is => 'String',
    	is_optional => 1,
    	doc => 'Output graph directory if any.',
    },
    plot_graph => {
    	is => 'Boolean',
    	is_optional => 1,
    	default => 0,
    	doc => 'Whether to plot all the graphs anyway no matter it passes filters of U and T test.',
    },
    tumor_normal_only => {
    	is => 'Boolean',
    	is_optional => 1,
    	default => 0,
    	doc => 'Whether to compute only tumor-normal p value, or all of the three.',
    },
    select_for_graph => {
    	is => 'Boolean',
    	is_optional => 1,
    	default => 0,
    	doc => 'Whether to plot graph according to the U test and T test. If yes, need to give U-test-threshold and T-test-threshold, and output directory.',
    },
    T_test_threshold => {
    	is => 'String',
    	is_optional => 1,
    	default => 0.01,
    	doc => 'If select for graph is chosen, will do T test between tumor and normal in ROI. Only p value smaller than this value may have the chance to have the graph plotted.',
    },
    mean_difference_threshold => {
    	is => 'String',
    	is_optional => 1,
    	default => 0.8,
    	doc => 'If select for graph is chosen, it will check if the difference between the mean of tumor and normal is large than this threshold if it already passed the T test step. If yes, the graph will be plotted.',
    },    
    name => {
    	is => 'String',
    	is_optional => 1,
    	doc => 'The name of the data if not specified in the file.',
    }
    ]
};

sub help_brief {
    "Do paired U test given the position of region of interest and bam file. It can call graph module if plot-graph is 1 or it passes both U test and further T test."
}

sub help_detail {
    "This script will do T test to compute P values for the region of interest (ROI) and the flanking region for both tumor and normal. It will also do T test to compare the ROI between tumor and normal. If select-for-graph is true, the script will statistics if it passes the T test filter. It will call graph module to draw a graph for each data if plot-graph is 1 or it passes two filters (T test, mean). graph-directory must be given if there is any possible graph output."
}  
	
sub execute {
    my $self = shift;

    # process input arguments
    my $outputFile = $self->output_file;
    my $inputFile = $self->input_file;
    
    my $isPlotGraph = $self->plot_graph;
    my $outputFigDir = $self->graph_directory;
    my $isTumorNormal = $self->tumor_normal_only;
    
    # choose big enough region as the flanking one for u test
    my $slide = 1000;
    my $point_number = 30;
    my $interval = $point_number * $slide;
    
    # choosing graphs to plot in file format
    my $selectPlot = $self->select_for_graph;
    my $tTestThreshold = $self->T_test_threshold;   
    my $mean_diff_threshold = $self->mean_difference_threshold;
    
    my $name=$self->name;

	# whether to write to a temporary file (1: write to a temporary file)    
   	my $system_tmp = 1;
    
    # Process options.
    if($inputFile !~/\S+/ || $outputFile !~/\S+/){
	    die("Either the input or output file was not given, or the input file does not exist. Please type 'gmt copy-number u-test -h' to see the manual.\n");
	}
	if(($isPlotGraph == 1 || $selectPlot == 1) && $outputFigDir !~/\S+/){
		die("Plot graph is possible but no graph-directory is given. Please type 'gmt copy-number u-test -h' to see the manual.\n");
	}
    if($isPlotGraph == 1 || $selectPlot == 1){
	    `mkdir $outputFigDir` unless (-e "$outputFigDir");    
	}
	
    #test architecture to make sure bam-window program can run (req. 64-bit)
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }
    
    my $db = "ucsc";
	my $user = "mgg_admin";
	my $password = "c\@nc3r"; 
	my $dataBase = "DBI:mysql:$db:mysql2";
	my $dbh = DBI->connect($dataBase, $user, $password) ||
    die "ERROR: Could not connect to database: $! \n";

	# Global hashes to convert UCSC name to Hugo name
	my ( %UcscToHugo, %UcscToUniprot, %UniprotToHugo );   

    # Set above hashes
	setNameConversionHashes(\%UcscToHugo, \%UcscToUniprot, \%UniprotToHugo, $dbh); 

	my ( $geneTableQuery, $geneStatement);
	my ($ucsc_name, $chrStart, $chrStop);


	# Query.  Call with $chr, $start, $end
	$geneTableQuery = "SELECT name, txStart, txEnd
                   FROM knownGene
                   WHERE chrom = ? && txEnd >= ? && txStart <= ?
                   ORDER BY txStart";
	$geneStatement = $dbh->prepare($geneTableQuery) ||
	die "Could not prepare statement '$geneTableQuery': $DBI::errstr \n";

   	open FILE_out, ">$outputFile" or die $!;
   	print FILE_out "name,chromosome,start-position,end-position,tumor-normal-t-value,tumor-p-value,normal-p-value,mean_tumor,mean_normal,gene-hit,plot-graph\n";
   	close FILE_out;
   	open FILE_out, ">>$outputFile" or die $!;
   	open FILE, "<$inputFile" or die $!;
   	while (my $line = <FILE>) {
   		chomp $line;
    	my ($name, $bam_tumor, $bam_normal, $chr, $start, $end) = split(/\,/,$line); 
		if(($bam_tumor !~/\S+/ && $bam_normal !~/\S+/) || $chr !~/\S+/ || $start !~/\d+/ || $end !~/\d+/){
		    die("Input file given but either the chromosome or positions or bam files are not given. Please type 'gmt copy-number u-test -h' to see the manual.\n");
		}
    	print "$chr\t$start\t$end\t";
    	if($end - $start > 1000000){
    		print FILE_out "$name,$chr,$start,$end,NA,NA,NA,NA,NA,NA,NA\n";
    		next;
    	}
    	
    	# define parameters
    	my $p1 = "NA";
    	my $p2 = "NA";
    	my $p3 = "NA";
    	my $p4 = "NA";
    	my $mean_tumor = "NA";
    	my $mean_normal = "NA";
    	my $gene_hit = 0;
    	my $select = 0;
    	
    	# check if gene hit		    
	    my @hasGene;		
	    my %hugoNames;
		$geneStatement->execute("chr$chr", $start, $end) || die "Could not execute statement '$geneTableQuery' with ($chr, $start, $end): $DBI::errstr \n";
		while ( ($ucsc_name, $chrStart, $chrStop) = $geneStatement->fetchrow_array() ){
			push(@hasGene, $ucsc_name);
		}
    
    	my $hasGene_count = @hasGene;
    	print "$hasGene_count\n";
	    if ( $hasGene_count ) {
	    	$gene_hit = 1;
	    } 
	    else {
	    	print FILE_out "$name,$chr,$start,$end,$p3,$p1,$p2,$mean_tumor,$mean_normal,$gene_hit,$select\n";
	    	next;
	    }
	    
   	    # read the neighbor boundary
   	    my $interval_;
   	    if($end - $start < $interval){
   	    	$interval_ = $interval;
   	    }
   	    else{
   	    	$interval_ = $end - $start;
   	    }
	    my $neighbor1_left = $start - $interval_;
	    my $neighbor1_right = $start - 1;
	    my $neighbor2_left = $end + 1;
	    my $neighbor2_right = $end + $interval_;
	    my $length = get_length($bam_tumor, $bam_normal, $chr);	    
		if($neighbor1_left < 0){
    		$neighbor1_left = 0;
    	}
    	if($neighbor1_right < 0){
    		$neighbor1_right = 0;
    	}
	    if($neighbor2_left > $length){
	    	$neighbor2_left = $length;
	    }
	    if($neighbor2_right > $length){
	    	$neighbor2_right = $length;
	    }
#print "$neighbor1_left\t$neighbor1_right\t$neighbor2_left\t$neighbor2_right\t$start\t$end\t$interval_\n";
    	# prepare the files   
	    my $tmp_in_tumor = "NA";
	    my $tmp_outL_tumor = "NA";
	    my $tmp_outR_tumor = "NA";
	    my $tmp_in_normal = "NA";
	    my $tmp_outL_normal = "NA";
	    my $tmp_outR_normal = "NA";
	    my $tmp_in_tumor_name_normalized = "NA";
	    my $tmp_in_normal_name_normalized = "NA";
	    
	    if($system_tmp == 1){ 
	    	if($bam_tumor =~ /\S+/){       
    		    my $tmp_in_tumor_normalized = File::Temp->new();
   			    $tmp_in_tumor_name_normalized = $tmp_in_tumor_normalized -> filename;
    		    my $tmp_in_tumor_ = File::Temp->new();
   			    $tmp_in_tumor = $tmp_in_tumor_ -> filename;
    		    my $tmp_outL_tumor_ = File::Temp->new();   	    
   			    $tmp_outL_tumor = $tmp_outL_tumor_ -> filename;
    	   	    my $tmp_outR_tumor_ = File::Temp->new();   	    
   			    $tmp_outR_tumor = $tmp_outR_tumor_ -> filename;
   		    }
    		if($bam_normal =~ /\S+/){
	    	    my $tmp_in_normal_normalized = File::Temp->new();   	    
   			    $tmp_in_normal_name_normalized = $tmp_in_normal_normalized -> filename;
    	   	    my $tmp_in_normal_ = File::Temp->new();
   			    $tmp_in_normal = $tmp_in_normal_ -> filename;
    		    my $tmp_outL_normal_ = File::Temp->new();   	    
   			    $tmp_outL_normal = $tmp_outL_normal_ -> filename;
    	   	    my $tmp_outR_normal_ = File::Temp->new();   	    
   			    $tmp_outR_normal = $tmp_outR_normal_ -> filename;
   			}
   		}
   		else{
	    	if($bam_tumor =~ /\S+/){       
	   			$tmp_in_tumor_name_normalized = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_tumor.csv";
	   			$tmp_in_tumor = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_in.csv";
	    	    $tmp_outL_tumor = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outL.csv";
	    	    $tmp_outR_tumor = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outR.csv";   			
   		    }
    		if($bam_normal =~ /\S+/){
	   			$tmp_in_normal_name_normalized = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_normal.csv";   		
   				$tmp_in_normal = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_inN.csv";
    		    $tmp_outL_normal = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outLN.csv";
    		    $tmp_outR_normal = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outRN.csv";   			
    		}
   		}	
    	     	
    	if($bam_tumor =~ /\S+/){
	    	$p1 = t_test($bam_tumor, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_tumor, $tmp_outL_tumor, $tmp_outR_tumor, $tmp_in_tumor_name_normalized, 1-$isTumorNormal);
	    }
    	if($bam_normal =~ /\S+/){
	    	$p2 = t_test($bam_normal, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_normal, $tmp_outL_normal, $tmp_outR_normal, $tmp_in_normal_name_normalized, 1-$isTumorNormal);    	
	    }
    	if($bam_tumor =~/\S+/ && $bam_normal =~ /\S+/){
	    	$p3 = t_test_ROI($bam_tumor, $bam_normal, $tmp_in_tumor_name_normalized, $tmp_in_normal_name_normalized);
	    }
    	
    	# do filtering
    	if($bam_tumor =~ /\S+/){
			$mean_tumor = mean_ROI($tmp_in_tumor_name_normalized);
		}
		if($bam_normal =~ /\S+/){
			$mean_normal = mean_ROI($tmp_in_normal_name_normalized);
		}
    	if($bam_tumor =~ /\S+/ && $bam_normal =~ /\S+/ && $selectPlot == 1 && $p3 < $tTestThreshold){
			if(abs($mean_tumor - $mean_normal) > $mean_diff_threshold){
				$select = 1;
			}
    	}
    	
    	# if not go through because of R defect, do it again.
    	my $count = 0;
    	if($bam_tumor =~ /\S+/ && $bam_normal =~ /\S+/){
	    	while(($p3 !~/\S+/ || $p1 !~/\S+/ || $p2 !~/\S+/ || $mean_tumor !~/\S+/ || $mean_normal !~/\S+/) && $count < 3){
    			$select = 0;
    			$count ++;
    			if($count == 1){
    				print "$chr,$start,$end:Missing values, will do it again!\n";
    			} else {
    				print "$chr,$start,$end:Missing values again, will do it again!\n";
    			}
    		
		    	$p1 = t_test($bam_tumor, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_tumor, $tmp_outL_tumor, $tmp_outR_tumor, $tmp_in_tumor_name_normalized, 1-$isTumorNormal);
    			$p2 = t_test($bam_normal, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_normal, $tmp_outL_normal, $tmp_outR_normal, $tmp_in_normal_name_normalized, 1-$isTumorNormal);    	
    	
    			$p3 = t_test_ROI($bam_tumor, $bam_normal, $tmp_in_tumor_name_normalized, $tmp_in_normal_name_normalized);
    	
    			# do filtering
				$mean_tumor = mean_ROI($tmp_in_tumor_name_normalized);
				$mean_normal = mean_ROI($tmp_in_normal_name_normalized);
    			if($selectPlot == 1 && $p3 < $tTestThreshold){
					if(abs($mean_tumor - $mean_normal) > $mean_diff_threshold){
						$select = 1;
					}
    			}
    		}    	
    	}
    		
    	if($bam_tumor =~ /\S+/ && $bam_normal =~ /\S+/ && $p3 !~/\S+/ || $p1 !~/\S+/ || $p2 !~/\S+/ || $mean_tumor !~/\S+/ || $mean_normal !~/\S+/){
    		$select = 0;
    		print "$chr,$start,$end:Still missing values, will leave it like what is the best!\n";
    	}
    	
		# write to the output file
		print FILE_out "$name,$chr,$start,$end,$p3,$p1,$p2,$mean_tumor,$mean_normal,$gene_hit,$select\n";		    	
			
		# deal with the plot: call R file directly
		if($isPlotGraph == 1 || $select == 1){
		
			# define the name of the picture with the genes
   			my $tmp_name;
    		if($system_tmp == 1){
			    my $tmp_ = File::Temp->new();
    			$tmp_name = $tmp_ -> filename;
		    }
		    else{
    			$tmp_name = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_name.csv";
		    }
	
			my $picName;
			# deal with the name	    
	        if($name !~/\S+/){
    		    $picName = $outputFigDir . "/Chr" . $chr . "_" . $start;
    		}
		    else{
        		$picName = $outputFigDir . "/". $name . "_chr" . $chr . "_" . $start;
		    }
		      
			# add the hugo gene name	
			my $count_geneName = 0;
			while ($hasGene[$count_geneName]){
			    my $hugo = convertUcscToHugo(\%UcscToHugo, \%UcscToUniprot, \%UniprotToHugo, $hasGene[$count_geneName]);
			    if ( defined $hugo ) { $hugoNames{$hugo} = 1; } 
			    $count_geneName ++;
			}    
		    
			if ( %hugoNames && scalar(keys %hugoNames) >= 1 ) { 
			    foreach (keys %hugoNames) { $picName = $picName . "_" . $_; }
			} 
			
		    $picName = $picName .  ".png";
    		    
		    my $isTitle = 1;
		    my $isSubTitle = 1;
		    my $isAnnotation = 0;
		    my $seg_file="";
		    my $rep_file="";
		    my $dgv_file="";
		    my $gene_file="";
		    my $array="";
		    my $tmp_pileup_name="";
		    my $tmp_pileup_name_n="";
		    my $isArray = 0;
		    my $isSnp = 0;
		    my $isFixYAxisLimit = 1;
    
		    open FILE_name, ">", $tmp_name or die $!;
		    print FILE_name "$tmp_in_tumor\t$tmp_outL_tumor\t$tmp_outR_tumor\t$tmp_in_normal\t$tmp_outL_normal\t$tmp_outR_normal\n$name\t$picName\t$chr\t$isTitle\t$isSubTitle\t$isAnnotation\n$start\t$end\t$neighbor1_left\t$neighbor1_right\t$neighbor2_left\t$neighbor2_right\n$seg_file\t$rep_file\t$dgv_file\t$gene_file\t$array\t$isArray\n$isSnp\t$tmp_pileup_name\t$tmp_pileup_name_n\t$isFixYAxisLimit\t\t\n";
		    close FILE_name;

		    my $command = qq{readcount(name="$tmp_name")};
		    my $library = "CN_graph.R";
		    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
		    $call -> execute;
		}
    }
	close FILE;	    	    
	close FILE_out;	
   	$dbh->disconnect();
}

sub u_test{
	my ($bam, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_name, $tmp_outL_name, $tmp_outR_name, $tmp_in_name_normalized, $uTest) = @_;

	if($bam eq ""){
		return "NA";
	}
	
	# read count for the three anyway
    write_read_count($bam, $chr, $start, $end, $tmp_in_name, $slide);
    write_read_count($bam, $chr, $neighbor1_left, $neighbor1_right, $tmp_outL_name, $slide);  
    write_read_count($bam, $chr, $neighbor2_left, $neighbor2_right, $tmp_outR_name, $slide);

    my $library = "tests.R";

    # normalization
    my $command = qq{normalize(name1="$tmp_in_name",nameL="$tmp_outL_name",nameR="$tmp_outR_name",normalizedFile="$tmp_in_name_normalized")};
    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    $call -> execute;

    # p value if U test
    if($uTest == 1){
        # temporary u test file
	    my $tmp_outAll_name = "NA";
	    my $tmp_outAll = File::Temp->new();
	   	$tmp_outAll_name = $tmp_outAll -> filename;

		# u test
	    $command = qq{utest(name1="$tmp_in_name",name_outL="$tmp_outL_name",name_outR="$tmp_outR_name",nameAll="$tmp_outAll_name")};
    	$call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    	$call -> execute;
    	
    	# read u test
	    open FILE_tmp, "<$tmp_outAll_name" or die $!;
	    my $line = <FILE_tmp>;
	    close FILE_tmp;
	    return $line;        
	}
	else{
		return "NA";
	}
}
  
sub u_test_ROI{
	my ($bam_tumor, $bam_normal, $tmp_tumor_name, $tmp_normal_name) = @_;
	if($bam_tumor eq "" || $bam_normal eq ""){
		return "NA";
	}
	
    my $tmp_outAll_name = "NA";
    my $system_tmp = 1;    
    if($system_tmp == 1){        
   	    my $tmp_outAll = File::Temp->new();
   	    $tmp_outAll_name = $tmp_outAll -> filename;
    }
    else{
  	    $tmp_outAll_name = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outAll_u_test.csv";
    }
    
    open FILE_tmp, ">$tmp_outAll_name" or die $!;
    close FILE_tmp;
                       
    my $command = qq{utest_ROI(name1="$tmp_tumor_name",name2="$tmp_normal_name",nameAll="$tmp_outAll_name")};
    my $library = "tests.R";
    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    $call -> execute;
    open FILE_tmp, "<$tmp_outAll_name" or die $!;
    my $line = <FILE_tmp>;
    return $line;        	
}
   
sub t_test{
	my ($bam, $chr, $start, $end, $neighbor1_left, $neighbor1_right, $neighbor2_left, $neighbor2_right, $slide, $tmp_in_name, $tmp_outL_name, $tmp_outR_name, $tmp_in_name_normalized, $tTest) = @_;

	if($bam eq ""){
		return "NA";
	}
	
	# read count for the three anyway
    write_read_count($bam, $chr, $start, $end, $tmp_in_name, $slide);
    write_read_count($bam, $chr, $neighbor1_left, $neighbor1_right, $tmp_outL_name, $slide);  
    write_read_count($bam, $chr, $neighbor2_left, $neighbor2_right, $tmp_outR_name, $slide);

    my $library = "tests.R";

    # normalization
    my $command = qq{normalize(name1="$tmp_in_name",nameL="$tmp_outL_name",nameR="$tmp_outR_name",normalizedFile="$tmp_in_name_normalized")};
    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    $call -> execute;

    # p value if T test
    if($tTest == 1){
        # temporary t test file
	    my $tmp_outAll_name = "NA";
	    my $tmp_outAll = File::Temp->new();
	   	$tmp_outAll_name = $tmp_outAll -> filename;

		# t test
	    $command = qq{ttest(name1="$tmp_in_name",name_outL="$tmp_outL_name",name_outR="$tmp_outR_name",nameAll="$tmp_outAll_name")};
    	$call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    	$call -> execute;
    	
    	# read u test
	    open FILE_tmp, "<$tmp_outAll_name" or die $!;
	    my $line = <FILE_tmp>;
	    close FILE_tmp;
	    return $line;        
	}
	else{
		return "NA";
	}
}
   
sub t_test_ROI{
	my ($bam_tumor, $bam_normal, $tmp_tumor_name, $tmp_normal_name) = @_;
	if($bam_tumor eq "" || $bam_normal eq ""){
		return "NA";
	}
	
    my $tmp_outAll_name = "NA";
    my $system_tmp = 1;    
    if($system_tmp == 1){        
   	    my $tmp_outAll = File::Temp->new();
   	    $tmp_outAll_name = $tmp_outAll -> filename;
    }
    else{
  	    $tmp_outAll_name = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outAll_t_test.csv";
    }
    
    open FILE_tmp, ">$tmp_outAll_name" or die $!;
    close FILE_tmp;
                       
    my $command = qq{ttest_ROI(name1="$tmp_tumor_name",name2="$tmp_normal_name",nameAll="$tmp_outAll_name")};
    my $library = "tests.R";
    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    $call -> execute;
    open FILE_tmp, "<$tmp_outAll_name" or die $!;
    my $line = <FILE_tmp>;
    return $line;        	
}
    
sub mean_ROI{
	my ($tmp_name) = @_;
	
    my $tmp_outAll_name = "NA";
    my $system_tmp = 1;    
    if($system_tmp == 1){        
   	    my $tmp_outAll = File::Temp->new();
   	    $tmp_outAll_name = $tmp_outAll -> filename;
    }
    else{
  	    $tmp_outAll_name = "/gscuser/xfan/svn/perl_modules/Genome/Model/Tools/Xian/tmp_outAll_mean.csv";
    }
    
    open FILE_tmp, ">$tmp_outAll_name" or die $!;
    close FILE_tmp;
                       
    my $command = qq{mean_ROI(name="$tmp_name",nameAll="$tmp_outAll_name")};
    my $library = "tests.R";
    my $call = Genome::Model::Tools::R::CallR->create(command=>$command, library=>$library);
    $call -> execute;
    open FILE_tmp, "<$tmp_outAll_name" or die $!;
    my $line = <FILE_tmp>;
    return $line;        	
}

sub write_read_count {
    my ($bam, $chr, $start, $end, $tmp_in_name, $slide) = @_;
    my @command = `samtools view $bam $chr:$start-$end`;
    open FILE_readcount, ">", $tmp_in_name or die $!;
    close FILE_readcount;
    open FILE_readcount, ">>", $tmp_in_name or die $!;

    # write read count
    my $current_window = $start;
    my $readcount_num = 0;
    for(my $i = 0; $i < $#command; $i ++ ) {
        my $each_line = $command[$i];
        my ($tmp1, $tmp2, $chr_here, $pos_here, $end_here,) = split(/\t/, $each_line);
        if($pos_here > $current_window + $slide){
        	if($readcount_num != 0){
            	print FILE_readcount "$chr\t$current_window\t$readcount_num\n";
            }
            if($current_window + $slide > $end){
                last;
            }
#            $current_window += $slide;
            if($pos_here > $current_window + 2*$slide){
            	$readcount_num = 0;
            }
            else{
	            $readcount_num = 1;
	        }
			$current_window += $slide;	        
        }
        else{
            $readcount_num ++;
        }
    } 
    if($current_window + $slide < $end){
    	# since it has not the chance to write the previous one
    	print FILE_readcount "$chr\t$current_window\t$readcount_num\n";
    	while($current_window + $slide < $end){
    		$current_window += $slide;
    		print FILE_readcount "$chr\t$current_window\t0\n";
    	}
    }
    close FILE_readcount;
    return;
}
    
sub get_length{
	my ($bam_tumor, $bam_normal, $chr) = @_;
	
	# get the length of the chromosome
	my $length;
	my $length1;
	my $length2;
	if(-e "$bam_tumor"){
		$length1 = readChrLength($bam_tumor, $chr);
#		print "$length1";
	}	
	if(-e "$bam_normal"){
		$length2 = readChrLength($bam_normal, $chr);
	}
	if(-e "$bam_tumor" && -e "$bam_normal"){
		if($length1 < $length2){
			$length = $length1;
		}
		else{
			$length = $length2;
		}
	}
	elsif(-e "$bam_tumor"){
			$length = $length1;
	}
	else{
			$length = $length2;
	}
	
	return $length;
}    

sub readChrLength{
    my ($bam, $chr) = @_;
	my @command = `samtools view -H $bam`;

    my $length;

	# skip the first two lines (first line: EOF; second line: header of the file)
	for(my $i = 1; $i < $#command; $i ++ ) {
        my $each_line = $command[$i];
        my ($tmp1, $chromosome_col, $length_col) = split(/\t/, $each_line);
        my ($tmp2, $chromosome) = split(/:/, $chromosome_col);
        if($chromosome eq $chr){
        	chomp $length_col;
        	(my $tmp3, $length) = split(/:/, $length_col);
	        last;
	    }
	}
	return $length;
}

sub readTable {
    my ($dbh, $table, $myChr, $myStart, $myStop, $myAnoFile, $geneTableQuery) = @_;
    # query
    # my $geneTableQuery = "SELECT chrom, chromStart, chromEnd FROM $table";
    my $geneStatement = $dbh->prepare($geneTableQuery) || die "Could not prepare statement '$geneTableQuery': $DBI::errstr \n";

    # execute query
    my ($chr, $chrStart, $chrStop);
    
    my $subString = ",";

    open FILE, ">", $myAnoFile or die $!;
    print FILE "Start\tEnd\n";
    close FILE;
    open FILE, ">>", $myAnoFile or die $!;
    $geneStatement->execute() || die "Could not execute statement for table knownGene: $DBI::errstr \n";
    while ( ($chr, $chrStart, $chrStop) = $geneStatement->fetchrow_array() ) {
        if($chr eq "chr".$myChr && $chrStart <= $myStop && $chrStop >= $myStart){ # overlap
            if($chrStart < $myStart){
                #$chrStart = $myStart;
                my $iIndex = index($myStart, $subString);
                if($iIndex >= 1){
                	$chrStart = substr($myStart, 0, $iIndex-1);
                }
                else{
                	$chrStart = $myStart;
                }
            }
            if($chrStop > $myStop){
                $chrStop = $myStop;
            }
            print FILE "$chrStart\t$chrStop\n";
        }
    }
    close FILE;
}


sub convertUcscToHugo {
    # Input: UCSC gene (transcript) ID
    # Return: Hugo ID or undef
    # First try UCSC -> Hugo look up based on current Hugo web site
    # Then do the indirect UCSC -> Uniprot -> Hugo

	my ( $UcscToHugo, $UcscToUniprot, $UniprotToHugo, $ucscId) = @_;     

    if ( defined $UcscToHugo->{$ucscId} ) { return $UcscToHugo->{$ucscId}; }

    if ( defined $UcscToUniprot->{$ucscId} &&
	 defined $UniprotToHugo->{$UcscToUniprot->{$ucscId}} ) {
	
	return $UniprotToHugo->{$UcscToUniprot->{$ucscId}};
    }

    # Can not unambiguously convert UCSC id to Hugo
    return undef;
}




sub setNameConversionHashes {
    # Set look-up hash for UCSC gene name to Hugo gene name
    # and Uniprot ID to Hugo ID
    #   UCSC -> Hugo should be unambiguous and will use the current UCSC ids 
    # Many of the current UCSC ids are not the ones in the build 36 database from 2006
    # Uniprot -> Hugo has some one to many (e.g.  P62158 => CALM1 and P62158 => CALM2; these are on different chromosomes)
    # If the Uniprot -> Hugo is ambiguous, don't use it.

	my ( $UcscToHugo, $UcscToUniprot, $UniprotToHugo, $dbh) = @_; 
    my ( $hugoGeneRef, $uniprotId, $hugo, $ucscId, %uniprotAmbiguous );
    
    $hugoGeneRef = HugoGeneMethods::makeHugoGeneObjects();
    foreach  $hugo ( keys %{$hugoGeneRef} ) {
	$ucscId = $$hugoGeneRef{$hugo}->ucsc();
	if ( defined $ucscId && $ucscId ne "" ) {
	    if ( defined $UcscToHugo->{$ucscId} && $UcscToHugo->{$ucscId} ne $hugo ) {
		print "WARNING: UCSC $ucscId => $hugo and $ucscId => $UcscToHugo->{$ucscId} \n";
	    }
	    $UcscToHugo->{$ucscId} = $hugo;
	}

	$uniprotId = $$hugoGeneRef{$hugo}->uniprot();
	if ( defined $uniprotId && $uniprotId ne "" ) {
	    if ( defined $UniprotToHugo->{$uniprotId} && $UniprotToHugo->{$uniprotId} ne $hugo ) {
		# this is ambiguous.  Add this Hugo ID to the rest
		$uniprotAmbiguous{$uniprotId} = 1;
		$UniprotToHugo->{$uniprotId} .= "|$hugo";
	    } else {
		$UniprotToHugo->{$uniprotId} = $hugo;
	    }
	}    
    }

    # Remove ambiguous Uniprot -> Hugo
    foreach (keys %uniprotAmbiguous ) { delete $UniprotToHugo->{$_}; }
        
    # Make hash for UCSC gene name to Uniprot ID based on UCSC database query
    # Doing this allows us to use old, build 36 UCSC names (in contrast of the real-time query to HGNC)
    # Then get Hugo from %uniprotToHugo
    my $query = "SELECT kgID, displayID FROM kgProtAlias";
    my $statement = $dbh->prepare($query)  ||
	die "Could not prepare statement '$query': $DBI::errstr \n";
    $statement->execute() ||
	die "Could not execute statement '$query': $DBI::errstr \n";
    while (  ($ucscId, $uniprotId ) = $statement->fetchrow_array() ) {
	if ( defined $UcscToUniprot->{$ucscId} && $UcscToUniprot->{$ucscId} ne $uniprotId ) {
	    print "OOPS: ambiguous UCSC to uniprot -- $ucscId, $uniprotId, $UcscToUniprot->{$ucscId} \n";
	} else {
	    $UcscToUniprot->{$ucscId} = $uniprotId;
	}
    }
}

