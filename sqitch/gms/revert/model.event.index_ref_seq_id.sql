-- Revert model.event.index_ref_seq_id

BEGIN;

DROP INDEX model.event_ref_seq_index;

COMMIT;
