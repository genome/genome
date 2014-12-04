-- Deploy remove_model_groups_from_analysis_project
-- requires: config_analysis_project

BEGIN;

  ALTER TABLE config.analysis_project DROP COLUMN model_group_id;

COMMIT;
