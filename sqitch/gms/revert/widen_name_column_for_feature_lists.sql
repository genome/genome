-- Revert widen_name_column_for_feature_lists

BEGIN;

  ALTER TABLE model.feature_list ALTER COLUMN name TYPE varchar(200);

COMMIT;
