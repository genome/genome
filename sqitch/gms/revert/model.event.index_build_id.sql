-- Revert model.event.index_build_id

BEGIN;

DROP INDEX model.event_build_id_index;

COMMIT;
