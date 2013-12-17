-- Revert model.event.index_event_type_model_id

BEGIN;

DROP INDEX model.event_type_model_id_index;

COMMIT;
