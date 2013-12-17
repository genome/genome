-- Revert model.event.index_parent_event_id

BEGIN;

DROP INDEX model.event_parent_event_index;

COMMIT;
