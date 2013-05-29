package Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFilesWorker;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFilesWorker {
    is  => ['Command'],
    has => [
        alignment_info => {
            is  => 'String',
            is_input => 1,
            doc => 'The directory containing the Maq map files.',
        },
        working_directory => {
            is => 'String',
            is_input => 1,
            doc => 'The working directory of the tool.',
        },
        ref_list => {
            is => 'String',
            is_input => 1,
            doc => 'The ref_list used for SamToBam.',
            example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai'],
            is_optional => 1,
        },
    ],
    has_param => [
           lsf_resource => {
           default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=2000]',
           },
    ],

};

sub help_brief {
    'Convert Maq map files into the TCGA format.';
}

sub help_detail {
    return <<EOS
    Convert Maq map files into the TCGA format.
EOS
}

sub execute {
    my $self = shift;
    my $working_dir = $self->working_directory;
    
    my $alignment_string = $self->alignment_info; 

    my @alignment_info = split(/\|/,$alignment_string);
    my $idid = shift @alignment_info; 
    my $alignment_directory = shift @alignment_info; 

    my $pid = getppid();
    my $log_path = "$working_dir/logs/convert_aligned_maps_$idid"."_".$pid.".txt";
    my $log_fh = Genome::Sys->open_file_for_writing($log_path);
    unless($log_fh) {
       $self->error_message("For $idid, failed to open output filehandle for: " .  $log_path );
       die "For $idid, could not open file ".$log_path." for writing.";
    } 

    print $log_fh "Alignment string: $alignment_string\n";
    print $log_fh "Alignment instrument-data id: $idid\n";
    print $log_fh "Alignment dir: $alignment_directory\n";

    my $mapmerge_output_file = "$working_dir/maps/$idid.map"; 
    my $conversion_output_file = "$working_dir/aligned/$idid.sam"; 

    #if no map file exists
    if (!-s $mapmerge_output_file) {
        my @maps = <$alignment_directory/*.map>;

        if ( scalar(@maps) == 0 ) {
            print $log_fh "\nNo maps found in $alignment_directory.  Returning.";
            return 1; 
        }

        if ( scalar(@maps) == 1  ) {
            #if there is only one map, no merge needed
            my $single_map = shift(@maps);
            print $log_fh "\nFound $single_map. Not merging, just copying to $mapmerge_output_file";
            if ($single_map =~ m/.+all_sequences.map$/) {
                #in this case, the all_sequences.map exist, just copy it to the proper location for merging
                my $copy_rv = Genome::Sys->copy_file($single_map,$mapmerge_output_file);
                if ($copy_rv ne 1)  {
                    print $log_fh "\nFor $idid, error copying $single_map to $mapmerge_output_file";
                    die "For $idid, error copying $single_map to $mapmerge_output_file";
                }        
            } else {
                #there is only one map and it is not all_sequences.map
                print $log_fh "\nNo appropriate map files found at: $alignment_directory.  Quitting.";
                die "No appropriate map files found at $alignment_directory";
            }
        } else { 
            #merge all the maps 
            print $log_fh "\nJoining...\n". join("\n",@maps) ;
            my $mapmerge_tool = Genome::Model::Tools::Maq::Mapmerge->create(use_version=>'0.7.1',
                                                                        input_map_files=>\@maps,
                                                                        output_map_file=>$mapmerge_output_file);
        
            my $mm_rv = $mapmerge_tool->execute();
            print $log_fh "\nDone with merge.  Merge return value is $mm_rv";

            if ($mm_rv ne 1) {
                $self->error_message("For $idid, error during merge: $mm_rv");
                die "For $idid, error during merge: $mm_rv";
            }
        }

    } else {
        print $log_fh "\n$mapmerge_output_file already exists.  Skipping generation of this file.\n";
    } 

    if (!-s $conversion_output_file) {

       # my $convert_cmd = $ENV{GENOME_SW} . "/samtools/samtools-0.1.6/misc/maq2sam-long $mapmerge_output_file $idid > $conversion_output_file"; 
       my $map_to_bam = Genome::Model::Tools::Maq::MapToBam->create(
                    map_file    => $mapmerge_output_file,
                    lib_tag     => $idid,
                    ref_list    => $self->ref_list,
                    index_bam   => 0,
                    sam_only    => 1,
       );

        print $log_fh "\nExecuting map to sam conversion: ";
        print $log_fh "\ninput file: $mapmerge_output_file";
        print $log_fh "\nlib tag: $idid";

        #my $rv = Genome::Sys->shellcmd( cmd=>$convert_cmd, input_files=>[$mapmerge_output_file], output_files=>[$conversion_output_file] );
        print $log_fh "\nResult from map2sam conversion: $map_to_bam";
        if ( $map_to_bam ne 1 ) {
            print $log_fh "\nError from map2sam conversion.";
            return;
        } else {
            print $log_fh "\nMap2sam conversion succeeded.";
        } 

        my $mapmerge_output_file =~ s/\.map$/\.sam/;
        #mv the file to the correct place
        my $mv_cmd = "mv $mapmerge_output_file $conversion_output_file"; 
        print $log_fh "\nRunning mv command: $mv_cmd";
        my $mv_rv = Genome::Sys->shellcmd( cmd=>$mv_cmd, input_files=>[$mapmerge_output_file], output_files=>[$conversion_output_file] );

        if ( $mv_rv ne 1 ) {
            print $log_fh "\nError while moving $mapmerge_output_file to $conversion_output_file"; 
            return;
        } else {
            print $log_fh "\nMove succeeded.";
        }
 
    } else {
        print $log_fh "\n$conversion_output_file already exists.  Skipping generation of this file.\n";
    } 

    $log_fh->close;
    return 1; 

}
1;
