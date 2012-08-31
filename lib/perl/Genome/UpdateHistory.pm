package Genome::UpdateHistory; 


class Genome::UpdateHistory {
    is => ['Genome::Notable'],
    table_name => 'MG.GENOME_UPDATE_HISTORY',
    id_by => [
         subject_table_name => { is => 'VARCHAR2', len => 255, column_name => 'SUBJECT_TABLE_NAME'},
                 subject_id => { is => 'VARCHAR2', len => 255, column_name => 'SUBJECT_ID'},
        subject_column_name => { is => 'VARCHAR2', len => 255, column_name => 'SUBJECT_COLUMN_NAME'},
                  edit_date => { is => 'TIMESTAMP', column_name => 'EDIT_DATE'}
    ],
    has => [
        description => { is => 'VARCHAR2', len => 255, column_name => 'DESCRIPTION', doc => 'insert, update, or delete'},
        is_reconciled => {is => 'Boolean', default => 0},
    ],
    has_optional => [
        old_value => { is => 'VARCHAR2', len => 1000, column_name => 'OLD_VALUE'},
        new_value => { is => 'VARCHAR2', len => 1000, column_name => 'NEW_VALUE'},
        app_user    => { is => 'VARCHAR2', len => 255, column_name => 'APP_USER'},
        app_name    => { is => 'VARCHAR2', len => 255, column_name => 'APP_NAME'},
        oracle_user => { is => 'VARCHAR2', len => 255, column_name => 'ORACLE_USER'},
        oracle_session_id => { is => 'VARCHAR2', len => 15, column_name => 'ORACLE_SESSION_ID'},
    ],
    data_source => 'Genome::DataSource::GMSchema',
};




1;


