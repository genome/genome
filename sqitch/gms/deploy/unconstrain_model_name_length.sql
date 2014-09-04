-- Deploy unconstrain_model_name_length
-- requires: model_model

BEGIN;

  ALTER TABLE model.model ALTER COLUMN name TYPE text;

COMMIT;
