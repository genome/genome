-- Revert model.event.index_run_id

BEGIN;

DROP INDEX model.event_run_id_index;

COMMIT;
