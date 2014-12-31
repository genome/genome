-- Deploy drop_featurelist_file_id
-- requires: model_feature_list

BEGIN;

  ALTER TABLE model.feature_list DROP COLUMN file_id;

COMMIT;
