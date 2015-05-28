package Genome::InstrumentData::AlignmentResult::ReliesOnBwa;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::AlignmentResult::ReliesOnBwa {
    is_abstract => 1,
};

# Object relying on bwa should just find a corresponding Bwa index and symlink it. This is the
# best we can do within the existing framework if we don't want to recreate an
# identical index already created by the 'regular' bwa module.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;

    my $bwa_version = $class->bwa_version($refindex);

    $class->debug_message("$class version $bwa_version is looking for a bwa version $bwa_version index.");

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );

    my $bwa_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        users              => $refindex->_user_data_for_nested_results,
        reference_build_id => $refindex->reference_build_id,
        aligner_name       => 'bwa',
        #aligner_params    => $refindex->aligner_params, # none of the aligner params should affect the index step so I think this okay
        aligner_version    => $bwa_version,
        test_name          => Genome::Config::get('aligner_index_test_name'),
    );

    for my $filepath (glob($bwa_index->output_dir . "/*")){
        my $filename = File::Basename::fileparse($filepath);
        next if $filename eq 'all_sequences.fa';
        next if $filename eq 'all_sequences.dict';
        Genome::Sys->create_symlink($filepath, $staging_dir . "/$filename");
    }

    $bwa_index->add_user(
        label => 'uses',
        user  => $refindex
    );

    return 1;
}


1;
