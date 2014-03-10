-- Verify model.event.index_event_status_event_type_model_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_e_es_et_mi';

ROLLBACK;
