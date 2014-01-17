-- Deploy add_clia_fields_to_analysis_projects
-- requires: config_analysis_project

BEGIN;

  ALTER TABLE config.analysis_project ADD COLUMN is_cle boolean;
  ALTER TABLE config.analysis_project ADD COLUMN run_as varchar(64);

COMMIT;
