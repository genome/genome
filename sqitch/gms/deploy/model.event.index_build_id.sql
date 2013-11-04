-- Deploy model.event.build_id
-- requires: model_event

BEGIN;

CREATE INDEX event_build_id_index on model.event using btree (build_id);

COMMIT;
