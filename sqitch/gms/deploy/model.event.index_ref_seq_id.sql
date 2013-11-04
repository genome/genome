-- Deploy model.event.ref_seq_id
-- requires: model_event

BEGIN;

CREATE INDEX event_ref_seq_index on model.event using btree (ref_seq_id);

COMMIT;
