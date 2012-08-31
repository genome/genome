package Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Dedup;

use strict;
use warnings;

use Genome;

use File::Basename;
use IO::File;

class Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Dedup {
    is => ['Command'],
    has_input => [
        accumulated_alignments_dir => {
            is  => 'String',
            doc => 'Accumulated alignments directory.' 
        },
        library_alignments => {
            is  => 'String',
            doc => 'Hash of library names and related alignment files.' 
        },
        aligner_version => {
            is  => 'Text',
            doc => 'The maq read aligner version used',
        },
        dedup_verison  => {
            is  => 'Text',
            doc => 'The duplication handler version used',
        },
        dedup_params  => {
            is  => 'Text',
            doc => 'The duplication handler params used',
        },
        ref_list    => {
            is  => 'Text',
            doc => 'The ref list that MapToBam uses to convert to bam',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=2000]',
        }
    ],
    has_output => [
        output_file => { 
            is => 'String', 
            is_optional => 1, 
        },
        library_name => { 
            is => 'String', 
            is_optional => 1, 
        }
    ],
};


sub make_real_rmdupped_map_file {
    my ($self, $maplist, $library, $log_fh, $working_directory) = @_;

    print $log_fh "Library: ".$library." Maplist: ".$maplist."\n";

    my $final_file = $self->accumulated_alignments_dir .  "/" .  $library.".map";
    print $log_fh "Final file: ". $final_file."\n";

    if (-s "$final_file") {
        print $log_fh "Rmdup'd file exists: ".$final_file."\n";
    } 
    else {
        my $tmp_file = Genome::Sys->create_temp_file_path($library.".map" );
        print $log_fh "Rmdup'd file DOES NOT exist: ".$final_file."\n";
        my $aligner_version = $self->aligner_version;
        my $maq_cmd = "gmt maq vmerge --version=$aligner_version --maplist $maplist --pipe $tmp_file &";
        print $log_fh "Executing:  $maq_cmd"."\n";
        #system "$maq_cmd";
        Genome::Sys->shellcmd( cmd=>$maq_cmd );
        my $start_time = time;
        until (-p "$tmp_file" or ( (time - $start_time) > 100) )  {
            sleep(5);
        }
        unless (-p "$tmp_file") {
            $log_fh->close;
            die "Failed to make intermediate file for (library) maps $!";
        }
        print $log_fh "Streaming into file $tmp_file."."\n";
        
        my $working_file = $working_directory . '/' . $library . ".map";
        
        my $maq_pathname = Genome::Model::Tools::Maq->path_for_maq_version($self->dedup_version);
        my $maq_tool = $maq_pathname . ' rmdup';
        my $dedup_params = $self->dedup_params;
        $maq_tool .= " $dedup_params" if $dedup_params;

        my $cmd = $maq_tool. " " . $working_file . " " . $tmp_file;
        print $log_fh "Running $cmd"."\n";
        #my $rv = system($cmd);
        my $rv = Genome::Sys->shellcmd(cmd=>$cmd);
        if($rv != 1) {
            print $log_fh "Problem with maq rmdup: $!"."\n";
            $log_fh->close;
            return;
        }
        
        rename($working_file, $final_file);
    }	
    return $final_file;
}

sub execute {
    my $self = shift;
    my $pid  = getppid(); 
    
    #Use Path::Class::Dir to correctly handle relative path when accumulated_alignments_dir is a symlink
    my $log_dir = Path::Class::Dir->new($self->accumulated_alignments_dir)->parent->subdir('logs')->stringify;
    unless (-e $log_dir ) {
	    unless( Genome::Sys->create_directory($log_dir) ) {
            $self->error_message("Failed to create log directory for dedup process: $log_dir");
            return;
	    }
    }
 
    my $log_file = $log_dir.'/parallel_dedup_'.$pid.'.log';
    my $log_fh = Genome::Sys->open_file_for_writing($log_file);
    unless($log_fh) {
       $self->error_message("Failed to open output filehandle for: " .  $log_file );
       die "Could not open file ".$log_file." for writing: " .  Genome::Sys->error_message;
    } 

    my $now = UR::Time->now;
    print $log_fh "Executing Dedup.pm at $now"."\n";

    my @list;
    if ( ref($self->library_alignments) ne 'ARRAY' ) {
        push @list, $self->library_alignments; 		
    } 
    else {
        @list = @{$self->library_alignments};   	#the parallelized code will only receive a list of one item. 
    }

    my $working_directory = File::Temp->newdir( 
        "parallel_dedup_$pid-working-XXXXX",
        DIR     => $self->accumulated_alignments_dir, 
        CLEANUP => 1
    );

    # fix permissions on this temp dir so others can clean it up later if need be    
    chmod(0775,$working_directory);

    print $log_fh "Input library list length: ".scalar(@list)."\n";
    for my $list_item ( @list ) {
        my %hash = %{$list_item};    		#there will only be one name-value-pair in the hash: $library name -> @list of alignment file paths (maps)
        for my $library ( keys %hash ) {
            $self->library_name($library);
            my @library_maps = @{$hash{$library}};
            print $log_fh "key:>$library<  /  value:>".scalar(@library_maps)."<"."\n";

            my $final_library_maplist = $self->accumulated_alignments_dir . '/' . $library . '.maplist';
            my $working_library_maplist = $working_directory .'/' . $library . '.maplist';
            print $log_fh "Library Maplist File:" .$final_library_maplist."\n";
            #my $fh = IO::File->new($library_maplist,'w');
            
            my $maplist_fh = Genome::Sys->open_file_for_writing($working_library_maplist);
            unless ($maplist_fh) {
                print $log_fh "Failed to create filehandle for '$working_library_maplist': " . Genome::Sys->error_message . "\n";
                $log_fh->close; 
                return;
            }
            my $cnt=0;
            for my $input_alignment (@library_maps) {
                unless(-f $input_alignment) {
                    print $log_fh "Expected $input_alignment not found.  Quitting."."\n";
                    $log_fh->close;
                    $maplist_fh->close; 
                    return;
                }
                $cnt++;
                print $maplist_fh $input_alignment ."\n";
            }
            print $log_fh "Library $library has $cnt map files"."\n";
            $maplist_fh->close;

            rename($working_library_maplist, $final_library_maplist);

            # db disconnect prior to map merge
            if (Genome::DataSource::GMSchema->has_default_handle) {
                $self->status_message("Disconnecting GMSchema default handle.");
                Genome::DataSource::GMSchema->disconnect_default_dbh();
            }

            $now = UR::Time->now;
            print $log_fh ">>> Starting make_real_rmdupped_map_file() at $now for library: $library ."."\n";
            my $map_file =  $self->make_real_rmdupped_map_file($final_library_maplist, $library, $log_fh, $working_directory);
            $now = UR::Time->now;
            print $log_fh "<<< Completed make_real_rmdupped_map_file() at $now for library: $library ."."\n";

            unless($map_file) {
                print $log_fh "Something went wrong with 'make_real_rmdupped_map_file'"."\n";
                $log_fh->close;
                return;
            }

            ###############
            #Beginning Map-2-Bam conversion

            $now = UR::Time->now;
            print $log_fh ">>> Beginning MapToBam conversion at $now for library: $library ."."\n";
            print $log_fh "MapToBam inputs for library: $library"."\n";
            print $log_fh "maq_version: ".$self->aligner_version."\n";
            print $log_fh "map_file: ".$map_file."\n";
            print $log_fh "lib_tag: ".$library."\n"; 
            
            my $map_to_bam = Genome::Model::Tools::Maq::MapToBam->create(
                use_version => $self->dedup_version,
                map_file    => $map_file,
                lib_tag     => $library,
                ref_list    => $self->ref_list,
                fix_mate    => 0,
            );
            my $map_to_bam_rv = $map_to_bam->execute;
            unless ($map_to_bam_rv == 1) {
                print $log_fh "MapToBam failed for library: $library with return value: $map_to_bam_rv"."\n";
                $log_fh->close;
                return;
            }
            $now = UR::Time->now;
            print $log_fh "<<< Ending MapToBam conversion at $now for library: $library ."."\n";
    
        }#end library loop 
        
        print $log_fh "*** Dedup process completed ***";
    }#end parallelized item loop

    $log_fh->close;
    return 1;
} #end execute


1;
