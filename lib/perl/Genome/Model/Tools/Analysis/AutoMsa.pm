package Genome::Model::Tools::Analysis::AutoMsa;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Analysis::AutoMsa {
    is => 'Command',                       
    has => [ 
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human; option under construction",
	    is_optional  => 1,
	    default => 'human',
	},
	ensemblversion => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2; option under construction",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},
	ace_fof => {
            type  =>  'String',
            doc  => "optional; provide an fof of ace files to run auto analysis including the full path",
 	    is_optional  => 1,
        },
	project_dir => {
            type  =>  'String',
            doc  => "manditory; full path of project-dir",
  	    is_optional  => 1,
	},
	ace_file => {
            type  =>  'String',
            doc  => "manditory; ace file name",
  	    is_optional  => 1,
	},
	refseq_id => {
            type  =>  'String',
            doc  => "optional; default attempts to find this name in the acefile; name of the refseq in the ace file",
  	    is_optional  => 1,
        },
	refseq_fasta => {
            type  =>  'String',
            doc  => "optional; default attempts to find this file based on the refseq_id; refseq_id.c1.refseq.fasta file",
  	    is_optional  => 1,
        },
	cds => {
            type  =>  'String',
            doc  => "obsolete in v_6 ;otherwise optional; default attempts to find this file based on the refseq name; cds.gff file",
   	    is_optional  => 1,
	},
	domain => {
            type  =>  'String',
            doc  => "obsolete in v_6 ;otherwise optional; default attempts to find this file based on the refseq name; domain.gff file",
   	    is_optional  => 1,
        },
	snp => {
            type  =>  'String',
            doc  => "optional; default attempts to find this file based on the refseq name; snp.gff file",
   	    is_optional  => 1,
        },
	amplicon => {
            type  =>  'String',
            doc  => 'obsolete in v_6 ;otherwise optional; default attempts to find this file based on the refseq name; amplicon.gff file',
   	    is_optional  => 1,
        },
	poly_source_1 => {
            type  =>  'String',
            doc  => 'optional; default is 1; poly-source-1 sets start of grouping of traces within a sample to be evaluated',
   	    is_optional  => 1,
	    default => '1',
        },
	poly_source_2 => {
            type  =>  'String',
            doc  => 'optional; default is 20; poly-source-2 sets stop of grouping of traces within a sample to be evaluated',
    	    is_optional  => 1,
 	    default => '20',
	},
	poly_indel_source_1 => {
            type  =>  'String',
            doc  => 'optional; default is 1; poly-indel-source-1 sets start of grouping of traces within a sample to be evaluated',
   	    is_optional  => 1,
 	    default => '1',
	},
	poly_indel_source_2 => {
            type  =>  'String',
            doc  => 'optional; default is 2; poly-indel-source-2 sets stop of grouping of traces within a sample to be evaluated',
    	    is_optional  => 1,
	    default => '2',
	},
	indelgroup => {
            type  =>  'String',
            doc  => 'indelgroup sets het indels within n bp will be jointly analyzed by indelBayes as used with polyscan; Default set at 50',
   	    is_optional  => 1,
 	    default => '50',
	},
	pretty_source_1 => {
            type  =>  'String',
            doc  => 'optional; default is 1; pretty-source-1 defines the start boundary of a sample name',
    	    is_optional  => 1,
 	    default => '1',
	},
	pretty_source_2 => {
            type  =>  'String',
            doc  => 'optional; default is 20; pretty-source-2 defines the end boundary of a sample name',
    	    is_optional  => 1,
 	    default => '20',
	},
	no_run_sift => {
            type  =>  'String',
            doc  => 'optional; default runs runsift; use no-run-sift option and runsift won\'t run',
   	    is_optional  => 1,
        },
	no_snp_detector => {
	    type  =>  'Boolean',
            doc  => 'use -no-snp-detector if you do not want to run snp-detector as a part of this analysis, use -no-snp-detector',
   	    is_optional  => 1,
        },
	mail_me => {
	    type  =>  'Boolean',
            doc  => 'optional; default no mail; use mail-me option to get mailed when the analysis is complete',
    	    is_optional  => 1,
	},
	analysis_version => {
	    type  =>  'String',
	    doc  => 'optional; default runs msa_autoanalysis_v6.0; however, if you stipulate 5 with this option you will run /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v5.0test.pl or 4 to run /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v4.1test.pl',
   	    is_optional  => 1,
  	    default => '6',
	},
	no_force_genotyping => {
	    type  =>  'Boolean',
            doc  => "use no-force-genotyping option default action is to force genotype",  #tag-mps
    	    is_optional  => 1,
	},
	force_genotype_coords_file => {
            type  =>  'String',
            doc  => "optional provide a file of coordinates to force genotype; File format must be first column chromosome second column coordinate columns should be seperated by space or tab this option only works with the default msa_autoanalysis_v6.0 version",
   	    is_optional  => 1,
        },
	run_snp_detector => {
	    type  =>  'Boolean',
	    doc  => "use --run-snp-detector if you dont want to run snp-detector as a part of this analysis default msa_autoanalysis_v6.0 version is not to run snp-detector",
	    is_optional  => 1,
        },
	lsf_memory_requirement => {
            type  =>  'String',
            doc  => "optional provide gigs of resource to reserve for your lsf job as a number from 4 to 8; Default is 4",
	    is_optional  => 1,
  	    default => '4',
	},
	], 
};


sub help_brief {
    "Analyze reads in a Consed Ace file"
}

sub help_synopsis {
    return <<EOS
	gmt analysis auto-msa -h
EOS
}
sub help_detail { 
    return <<EOS 

This tool was designed to set up and launch auto analysis of a consed ace file on the blades.
It is manditory to use either an ace fof or an ace file and project dir. 

        See gmt analysis auto-msa --help for a full list of options

   minimun arguments for running the tool are -ace-file and -project-dir or -ace-fof 
EOS
}


sub execute {

    my $self = shift;

    my $ace_fof = $self->ace_fof;
    my $ace_file = $self->ace_file;
    my $project_dir = $self->project_dir;

    my $mail_me = $self->mail_me;

    my $send_it;
    my ($handle) = getpwuid($<);

    unless ($ace_fof || ($ace_file && $project_dir)) {
	$self->error_message( "\n\nIt is manditory to use either an ace_fof or and ace_file and project_dir."
			    . " invalid command ."
			    . " See gmt analysis auto-msa --help for options.\n\n");
	
	if ($mail_me) {$send_it =  qq(echo "It is manditory to use either an ace fof or an ace file and project dir see gmt analysis auto-msa --help for options. You'll need to restart your analysis." | mailx -s 'Invalid command. Anaysis did not start.' $handle);`$send_it`}
	return;
    }
    
    my $launch_dir=`pwd`;
    chomp($launch_dir);

    if ($ace_fof) {
	$send_it =  qq(echo "Could not open the ace_fof $ace_fof.  Check your ace_fof and restart your analysis." | mailx -s 'Invalid command. Anaysis did not start.' $handle);
	unless ($mail_me) {$send_it = '';}
	unless (-f $ace_fof) { print "\n\nCould not open the ace_fof $ace_fof. Check your ace_fof and restart your analysis.\n\n";`$send_it`;return; }
	open(FOF,$ace_fof);
	my $n = 0;
	my $r_n = 0;

	my @failed_to_launch;
	while (<FOF>) {
	    chomp;
	    my $ace = $_;
	    $n++;
	    my ($project_dir,$ace_file) = $ace =~ /^(\S+)\/edit\_dir\/(\S+)$/;
	    
	    chdir $launch_dir;
	    my $run = &set_up_and_run($project_dir,$ace_file,$self,$handle);
	    if ($run) {
		$r_n++;
	    } else {
		push (@failed_to_launch,$project_dir);
	    }
	}
	close FOF;
	if ($r_n == $n) {
	    print "\n\nAn analysis for each of the projects ($r_n) in your ace fof was submitted to the blades for processing.\n\n";
	    $send_it =  qq(echo "An analysis for each of the $r_n projects in your ace fof $ace_fof was submitted to the blades for processing.\n" | mailx -s 'All $n projects in $ace_fof were submitted to the blades for processing.' $handle);
	    if ($mail_me) {`$send_it`;}
	    return 1;
	} else {
	    my $failed_to_launch_projects = join '\t\n' , @failed_to_launch;
	    $self->error_message( "\n\nOnly $r_n of the $n projects in your ace fof were submitted to the blades for processing.  The following projects will need to be checked and restarted;\n\t$failed_to_launch_projects.\n\n"); 
	    $send_it =  qq(echo "Only $r_n of the $n projects in your ace fof was submitted to the blades for processing. The following projects will need to be checked and restarted;\n\t$failed_to_launch_projects." | mailx -s 'Only $r_n of the $n projects in $ace_fof were submitted to the blades for processing.' $handle);
	    if ($mail_me) {`$send_it`;}
	    return;
	}
    } else {
	my $run = &set_up_and_run($project_dir,$ace_file,$self,$handle);
	if ($run) {
	    print "\n\nYour analysis of $project_dir $ace_file was submitted to the blades for processing.\n\n";
	    $send_it =  qq(echo "Your analysis of $project_dir $ace_file was submitted to the blades for processing.\n" | mailx -s 'Your analysis of $project_dir $ace_file was submitted to the blades for processing.' $handle);
	    if ($mail_me) {`$send_it`;}
	    return 1;
	} else {
	    return;
	}
    }
}


sub set_up_and_run {

    my ($project_dir,$ace_file,$self,$handle) = @_;
    my $check = &check_ace_dir($ace_file,$project_dir,$handle,$self);
    unless ($check) {return;}

######################################################################################
###This section uses the ace object to find the refseq id
######################################################################################
    my $refseq_id = $self->{refseq_id};
    my $ace_fof = $self->ace_fof;
    my $analysis_version = $self->analysis_version;
    my $send_it;
    my $ed = "$project_dir/edit_dir";

    my $cds = $self->cds;
    my $domain = $self->domain;
    my $snp = $self->snp;
    my $amplicon = $self->amplicon;
    my $poly_source_1 = $self->poly_source_1;
    my $poly_source_2 = $self->poly_source_2;
    my $poly_indel_source_1 = $self->poly_indel_source_1;
    my $poly_indel_source_2 = $self->poly_indel_source_2;
    my $indelgroup = $self->indelgroup;
    my $pretty_source_1 = $self->pretty_source_1;
    my $pretty_source_2 = $self->pretty_source_2;
    my $no_run_sift = $self->no_run_sift;
    my $no_snp_detector = $self->no_snp_detector;
    my $mail_me = $self->mail_me;
    my $tag_mps = $self->no_force_genotyping;
    my $force_genotype_coords_file = $self->force_genotype_coords_file;
    my $run_snp_detector = $self->run_snp_detector;
    
    my $get_ref_id;
    if ($ace_fof) { $get_ref_id = 1; }
    unless ($refseq_id) { $get_ref_id = 1; }

    if ($get_ref_id) {
        my $Contig_number;

        my $ao = Genome::Model::Tools::Consed::AceReader->create(file => "$ed/$ace_file");
        die $self->error_message("Failed to create ace reader for file: $ed/$ace_file") if not $ao;

        while ( my $contig = $ao->next_contig ) {
            if (grep { /\.c1$/ } keys %{ $contig->{reads} }) {
                $Contig_number=$contig->{name};
            }
            foreach my $read_name (keys %{ $contig->{reads} }) {
                if ($read_name =~ /(\S+\.c1)$/) {
                    if ($refseq_id) {
                        $send_it =  qq(echo "More than one refseq_id was found. Please stipulate refseq_id as a command line option and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. More than one refseq_id was found.' $handle);
                        if ($mail_me) {`$send_it`;}
                        $self->error_message( "\n\nMore than one refseq_id was found. Please stipulate refseq_id as a command line option. See "
                            . App::Name->prog_name
                            . " --help for more command line parameters.\n\n" );
                        return;
                    } else {
                        $refseq_id=$read_name;
                    }
                }
            }
        }
    }

    unless ($refseq_id) {
        $send_it =  qq(echo "Could not identify the refseq_id. Please check the Ace file and stipulate refseq_id as a command line option and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not identify the refseq_id.' $handle);
        if ($mail_me) {`$send_it`;}
        print qq(\nCould not identify the refseq_id. Please check the Ace file and stipulate refseq_id as a command line option and restart your analysis of $ed/$ace_file.\n\n);
        return;
    }

    print "$refseq_id\n";

######################################################################################
######################################################################################

    #--- Check Current directory and verify the existance of files--- 

    my $cwd=`pwd`;
    chomp($cwd);
    
    if ($cwd ne $ed) {
	chdir $ed;
    }
    
    my $newcd = `pwd`; # eq $ed
    print "$newcd\n";
    
    my $refseq_fasta = $self->refseq_fasta;
    
    if ($refseq_fasta) {
	unless ($refseq_fasta && -e $refseq_fasta) {
	    $send_it =  qq(echo "Could not find the refseq fasta file provided $refseq_fasta. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the refseq fasta.' $handle);
	    if ($mail_me) {`$send_it`;}
	    $self->error_message( "\n\nCould not find the refseq fasta file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
	    return;
	}
	
    } else {
	$refseq_fasta = "$refseq_id.refseq.fasta";
	unless ($refseq_fasta && -e $refseq_fasta) {
	    $send_it =  qq(echo "Could not find the refseq fasta file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the refseq fasta.' $handle);
	    if ($mail_me) {`$send_it`;}
	    $self->error_message( "\n\nCould not find the refseq fasta file in the edit_dir. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
	    return;
	}
    }
    
    my $root = (split(/\.c1/,$refseq_id))[0]; #nnd = number name date
    
    unless ($analysis_version == 6) {
	my ($gene_id, $gene, $d_date) = (split(/\_/,$root))[0,1,2];
	
	if ($cds) {
	    unless ($cds && -e $cds) {
		$send_it =  qq(echo "Could not find the cds file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the cds file.' $handle);
		if ($mail_me) {`$send_it`;}
		$self->error_message( "\n\nCould not find the cds file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
		return;
	    }
	} else {
	    $cds = <$gene_id.$gene.*.$d_date.cds.gff>;
	    unless ($cds && -e $cds) {
		$cds = "$gene_id.$gene.*.$d_date.cds.gff";
		unless ($cds && -e $cds) {
		    my $tsp_root = (split(/\.c1/,$refseq_id))[0];
		    $cds = "$tsp_root.cds.gff";
		    ($gene_id, $gene, $d_date) = (split(/\-/,$tsp_root))[0,1,3]; ###redefined here but redefined incorrectly
		    #TSP_Project-0007157-Ensembl-36_35i
		    unless ($cds && -e $cds) {
			$send_it =  qq(echo "Could not find the cds file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the cds file.' $handle);
			if ($mail_me) {`$send_it`;}
			$self->error_message( "\n\nCould not find the cds file. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
			return;
			
		    }
		}
	    }
	}
	
	if ($cds) {
	    $root = (split(/\.cds.gff/,$cds))[0];
	    
	    if ($domain) {
		unless ($domain && -e $domain) {
		    $send_it =  qq(echo "Could not find the domain file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the domain file.' $handle);
		    if ($mail_me) {`$send_it`;}
		    $self->error_message( "\n\nCould not find the domain file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
		    return;
		}
	    } else {
		$domain = "$root.domain.gff";
		unless ($domain && -e $domain) {
		    $domain = <$gene_id.$gene.*.$d_date.domain.gff>;	
		    unless ($domain && -e $domain) {
			undef($domain);
			$self->error_message( "\n\nCould not find the domain file in the edit_dir. Autoanalysis will continue without it.\n\n"); 
		    }
		}
	    }
	    
	    if ($snp) {
		unless ($snp && -e $snp) {
		    $send_it =  qq(echo "Could not find the dbsnp file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the dbsnp file.' $handle);
		    if ($mail_me) {`$send_it`;}
		    $self->error_message( "\n\nCould not find snp file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
		    return;
		}
	    } else {
		$snp = "$root.snp.gff";
		unless ($snp && -e $snp) {
		    $snp = <$gene_id.$gene.*.$d_date.domain.gff>;
		    unless ($snp && -e $snp) {
			undef($snp);
			$self->error_message( "\n\nCould not find the dbsnp file in the edit_dir. Autoanalysis will continue without it.\n\n"); 
		    }
		}
	    }
	    
	    if ($amplicon) {
		unless ($amplicon && -e $amplicon) {
		    $send_it =  qq(echo "Could not find the amplicon file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the amplicon file.' $handle);
		    if ($mail_me) {`$send_it`;}
		    $self->error_message( "\n\nCould not find the amplicon file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n");
		    return;
		}
	    } else {
		$amplicon = "$root.amplicon.gff";
		unless ($amplicon && -e $amplicon) {
		    $amplicon = <$gene_id.$gene.*.$d_date.domain.gff>;
		    unless ($amplicon && -e $amplicon) {
			undef($amplicon);
			$self->error_message( "\n\nCould not find the amplicon file in the edit_dir. Autoanalysis will continue without it.\n\n"); 
		    }
		}
	    }
	}
    }


##Build the command line / running conditions
    
######################################################################
	    
    my @msa_auto_cmd;
    my @default;
    if ($refseq_fasta && -e $refseq_fasta) {
	
	unless ($analysis_version == 6) {
	    @msa_auto_cmd = ("--project-dir", "$project_dir", "--ace-file", "$ace_file", "--cds", "$cds", "--refseq-id", "$refseq_id", "--refseq-fasta", "$refseq_fasta", "--poly-source-1", "$poly_source_1", "--poly-source-2", "$poly_source_2", "--poly-indel-source-1", "$poly_indel_source_1", "--poly-indel-source-2", "$poly_indel_source_2", "--indelgroup", "$indelgroup", "--pretty-source-1", "$pretty_source_1", "--pretty-source-2", "$pretty_source_2");
	}
	
	@default  = ("--project-dir", "$project_dir", "--ace-file", "$ace_file", "--refseq-id", "$refseq_id", "--refseq-fasta", "$refseq_fasta", "--poly-source-1", "$poly_source_1", "--poly-source-2", "$poly_source_2", "--poly-indel-source-1", "$poly_indel_source_1", "--poly-indel-source-2", "$poly_indel_source_2", "--indelgroup", "$indelgroup", "--pretty-source-1", "$pretty_source_1", "--pretty-source-2", "$pretty_source_2");
	
	if ($domain) {
	    push (@msa_auto_cmd, "--domain");
	    push (@msa_auto_cmd, $domain);
	}
	if ($snp && -e $snp) {
	    push (@msa_auto_cmd, "--snp");
	    push (@msa_auto_cmd, $snp);
	    push (@default, "--dbsnp-file");
	    push (@default, $snp);
	    unless ($snp && -e $snp) { #added the extra check here for the benifit of v6
		$send_it =  qq(echo "Could not find the dbsnp file. Check your edit_dir and restart your analysis." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the dbsnp file.' $handle);
		if ($mail_me) {`$send_it`;}
		$self->error_message( "\n\nCould not find snp file provided. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
		return;}
	}
	if ($amplicon) {
	    push (@msa_auto_cmd, "--amplicon");
	    push (@msa_auto_cmd, $amplicon);
	}
	if ($no_run_sift) {
	    push (@msa_auto_cmd, "--no-run-sift");
	    push (@msa_auto_cmd, $no_run_sift);
	    push (@default, "--no-run-sift");
	}
	if ($no_snp_detector) {
	    push (@msa_auto_cmd, "--no-snp-detector");
	    push (@msa_auto_cmd, $no_snp_detector);
	}
	if ($mail_me) {
	    push (@msa_auto_cmd, "--mail-me");
	    push (@msa_auto_cmd, $mail_me);
	    push (@default, "--mail-me");
	}
	if ($tag_mps) {
	    push (@msa_auto_cmd, "--tag-mps");
	    push (@msa_auto_cmd, $tag_mps);
	    push (@default, "--tag-mps");
	}
	if ($force_genotype_coords_file && -e $force_genotype_coords_file) {
	    push (@default, "--force-genotype-coords");
	    push (@default, $force_genotype_coords_file);
	}
	if ($run_snp_detector) {
	    push (@default, "--run-snp-detector");
	}
	
	my $run_autoanalysis = join(' ', @msa_auto_cmd);
	my $run_v6 =  join(' ', @default);
	
######################################################################

	my $memory_requirement = $self->lsf_memory_requirement;
	unless ($memory_requirement =~ /^[4 5 6 7 8]$/) {
	    $self->error_message( "\n\nyour request for $memory_requirement as a memory requirement has been denied. Your job will run with 4G memory. See gmt analysis auto-msa -h for more details.\n");
	    $memory_requirement = 4;
	}
	my $lsf_memory_requirement = "-R 'select[mem>$memory_requirement" . "000] rusage[mem=$memory_requirement" . "000]' -M $memory_requirement" . "000000";
	    print "\n\nThese are the files and run parameter's your analysis will use.\n$lsf_memory_requirement $run_v6\n\n";

	if ($analysis_version == 6) {
	    print "\nYou will run the default version 6. This runs the deployed version msa_autoanalysis_v6.0.\n";
	    
	    system ("bsub -oo msaSTDout_V6 $lsf_memory_requirement msa_autoanalysis_v6.0 $run_v6");
	    print "\n\nThese are the files and run parameter's your analysis will use.\n$run_v6\n\n";
    
	} else {

	    if ($cds && -e $cds) {
		if ($analysis_version == 4) {
		    print "\nYou are running version 4 /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v4.1test.pl $run_autoanalysis\n";
		
		    system ("bsub -oo msaSTDout_v4.1 $lsf_memory_requirement /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v4.1test.pl $run_autoanalysis");
		    print "\nThe current test version is testing the lastest polyscan and combine_variants\nThese are the files and run parameter's your analysis will use.\n$run_autoanalysis\n\n";
		} elsif ($analysis_version == 5) {
		    print "\nYou are running version 5 /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v5.0test.pl $run_autoanalysis\n";
		    
		    system ("bsub -oo msaSTDout_v5 $lsf_memory_requirement /gscmnt/200/medseq/analysis/software/scripts/msa_autoanalysis_v5.0test.pl $run_autoanalysis");
		    print "\nThe current test version is testing the lastest polyscan and combine_variants\nThese are the files and run parameter's your analysis will use.\n$run_autoanalysis\n\n";

		}
	    } 
	}
	return 1;
    } else {
	$send_it =  qq(echo "Could not find the refseq fasta file. Check your edit_dir and restart your analysis of $ed/$ace_file." | mailx -s 'Anaysis did not start for $ed/$ace_file. Could not find the refseq fasta.' $handle);
	if ($mail_me) {`$send_it`;}
	$self->error_message( "\n\nCould not find the refseq fasta file in the edit_dir. Check your edit_dir and restart your analysis of $ed/$ace_file.\n\n"); 
	return;
    }
}

sub check_ace_dir {

    my ($ace_file,$project_dir,$handle,$self) = @_;
    my $mail_me = $self->mail_me;
    my $send_it =  qq(echo "You'll need to restart your analysis of $project_dir $ace_file ." | mailx -s 'Analysis of $ace_file failed.' $handle);
    my $ace_fof = $self->ace_fof;
    
    unless ($project_dir && -e "$project_dir") {
	$self->error_message( "\n\nCould not see the project directory.$project_dir. Your analysis of $project_dir $ace_file failed to start.\n\n");
	$send_it =  qq(echo "You'll need to restart your analysis of $project_dir $ace_file. Could not see the project directory." | mailx -s 'Analysis of $ace_file failed.' $handle);
	unless (qq($project_dir/edit_dir/$ace_file) && -e qq($project_dir/edit_dir/$ace_file)) { 
	    $self->error_message( "\n\nCould not find the ace file. $project_dir/edit_dir/$ace_file\n\n");
	    $send_it =  qq(echo "You'll need to restart your analysis of $project_dir $ace_file. Could not see the project directory and Could not find the ace file." | mailx -s 'Analysis of $ace_file failed.' $handle);
	}
	if ($mail_me) {`$send_it`;}
	return;
    }
    
    unless (qq($project_dir/edit_dir/$ace_file) && -e qq($project_dir/edit_dir/$ace_file)) { 
	$send_it =  qq(echo "You'll need to restart your analysis of $project_dir $ace_file. Could not find the ace file." | mailx -s 'Analysis of $ace_file failed.' $handle);
	if ($mail_me) {`$send_it`;}
	$self->error_message( "\n\nCould not find the ace file. Your analysis of $project_dir/edit_dir/$ace_file failed to start.\n\n");
	return;
    }
    return 1;
}
