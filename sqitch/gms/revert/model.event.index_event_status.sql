-- Revert model.event.index_event_status

BEGIN;

DROP INDEX model.event_status_index;

COMMIT;
