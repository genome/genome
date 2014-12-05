-- Deploy make_model_group_id_nullable
-- requires: config_analysis_project

BEGIN;

  ALTER TABLE config.analysis_project ALTER COLUMN model_group_id DROP NOT NULL;

COMMIT;
