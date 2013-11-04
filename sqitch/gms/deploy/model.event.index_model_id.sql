-- Deploy model.event.model_id
-- requires: model_event

BEGIN;

CREATE INDEX event_model_id_index on model.event using btree (model_id);

COMMIT;
