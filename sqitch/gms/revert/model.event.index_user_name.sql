-- Revert model.event.index_user_name

BEGIN;

DROP INDEX model.event_user_name_index;

COMMIT;
