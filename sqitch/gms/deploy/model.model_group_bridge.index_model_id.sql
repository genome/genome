-- Deploy model.model_group_bridge.model_id
-- requires: model_model_group_bridge

BEGIN;

CREATE INDEX model_group_bridge_model_id on model.model_group_bridge using btree (model_id);

COMMIT;
