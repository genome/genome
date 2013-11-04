-- Verify model.event.index_event_type_model_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_type_model_id_index';

ROLLBACK;
