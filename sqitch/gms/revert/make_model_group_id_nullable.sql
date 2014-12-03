-- Revert make_model_group_id_nullable

BEGIN;

  ALTER TABLE config.analysis_project ALTER COLUMN model_group_id SET NOT NULL;

COMMIT;
