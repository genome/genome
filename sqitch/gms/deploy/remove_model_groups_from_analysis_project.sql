-- Deploy remove_model_groups_from_analysis_project
-- requires: make_model_group_id_nullable

BEGIN;

  ALTER TABLE config.analysis_project DROP COLUMN model_group_id;

COMMIT;
