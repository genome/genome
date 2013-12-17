-- Deploy model.event.run_id
-- requires: model_event

BEGIN;

CREATE INDEX event_run_id_index on model.event using btree (run_id);

COMMIT;
