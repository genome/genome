-- Deploy model.event.event_type_model_id
-- requires: model_event

BEGIN;

CREATE INDEX event_type_model_id_index on model.event using btree (event_type, model_id);

COMMIT;
