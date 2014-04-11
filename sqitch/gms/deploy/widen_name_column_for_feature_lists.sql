-- Deploy widen_name_column_for_feature_lists
-- requires: model_feature_list

BEGIN;

  ALTER TABLE model.feature_list ALTER COLUMN name TYPE varchar(1500);

COMMIT;
