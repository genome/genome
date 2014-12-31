-- Revert drop_featurelist_file_id

BEGIN;

  ALTER TABLE model.feature_list ADD COLUMN file_id bigint;

COMMIT;
