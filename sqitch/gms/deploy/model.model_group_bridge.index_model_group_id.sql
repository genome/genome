-- Deploy model.model_group_bridge.model_group_id
-- requires: model_model_group_bridge

BEGIN;

CREATE INDEX model_group_bridge_model_group_id on model.model_group_bridge using btree (model_group_id);

COMMIT;
