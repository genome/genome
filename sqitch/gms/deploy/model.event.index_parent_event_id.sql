-- Deploy model.event.parent_event_id
-- requires: model_event

BEGIN;

CREATE INDEX event_parent_event_index on model.event using btree (parent_event_id);

COMMIT;
