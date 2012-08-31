package Genome::Library::Command::MergeDuplicates;
class Genome::Library::Command::MergeDuplicates {
    is => 'Genome::Command::Base',
    has_many => [
        libraries => { is => 'Genome::Library', shell_args_position => 1 },
    ],
};

sub execute {
    my $self = shift;
    my %name_to_libs;
    push @{$name_to_libs{$_->name}}, $_ for ($self->libraries);

    print "Found " . scalar(keys %name_to_libs) . " groups\n";

    while(my ($name, $libs) = each %name_to_libs){
        print scalar(@$libs) . " libraries have name $name\n";
        my $final_lib = shift @$libs;
        for my $old_lib (@$libs) {
            print "Moving data from library " . $old_lib->id . "\n";
            for my $data ($old_lib->instrument_data){
                print "Data " . $data->id . " moving to library " . $final_lib->id . "\n";
                $data->library($final_lib);
            }
            print "Deleting library " . $old_lib->id . "\n";
            $old_lib->delete;
        }
    }
    return 1;
}
