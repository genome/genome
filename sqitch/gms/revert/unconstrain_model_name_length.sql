-- Revert unconstrain_model_name_length

BEGIN;

  ALTER TABLE model.feature_list ALTER COLUMN name TYPE varchar(255);

COMMIT;
