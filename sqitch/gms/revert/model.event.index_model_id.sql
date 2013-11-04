-- Revert model.event.index_model_id

BEGIN;

DROP INDEX model.event_model_id_index;

COMMIT;
