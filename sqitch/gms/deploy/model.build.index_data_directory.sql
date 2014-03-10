-- Deploy model.build.data_directory
-- requires: model_build

BEGIN;

CREATE INDEX build_directory_index on model.build using btree (data_directory);

COMMIT;
