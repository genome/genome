-- Deploy model.feature_list.name
-- requires: model_feature_list

BEGIN;

CREATE INDEX feature_list_name_index on model.feature_list using btree (name);

COMMIT;
