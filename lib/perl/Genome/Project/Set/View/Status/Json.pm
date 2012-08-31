
package Genome::Project::Set::View::Status::Json;


class Genome::Project::Set::View::Status::Json {
    is => 'UR::Object::View::Default::Json',
    has => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'count',
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'json',
                    aspects => ['id', 'name', 'fixed_size_name','parts_count'],
                }
            ]
        }
    ]
};



1;



