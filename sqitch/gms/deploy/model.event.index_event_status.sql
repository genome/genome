-- Deploy model.event.event_status
-- requires: model_event

BEGIN;

CREATE INDEX event_status_index on model.event using btree (event_status);

COMMIT;
