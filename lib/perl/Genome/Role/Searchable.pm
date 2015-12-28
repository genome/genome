package Genome::Role::Searchable;

use Genome;
use UR::Role;

# A role indicating that this class creates a searchable item in the Solr system
role Genome::Role::Searchable { };

my $tried_loading_search_pm;
sub __import__ {
    my($role_name, $class_meta) = @_;
    my $class_name = $class_meta->class_name;

    # if the search engine is installed, configure its hooks
    unless ($tried_loading_search_pm++) {
        local($SIG{__WARN__});
        local($SIG{__DIE__});
        require Genome::Search;
    }

    # This ensures that the search system is updated when certain classes are updated
    # The search system is optional so it skips this if usage above fails
    if ($INC{"Genome/Search.pm"}) {
        Genome::Search->register_callbacks($class_name);
    }
}

1;
