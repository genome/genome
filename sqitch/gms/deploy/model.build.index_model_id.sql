-- Deploy model.build.model_id
-- requires: model_build

BEGIN;

CREATE INDEX build_model_index on model.build using btree (model_id);

COMMIT;
