-- Verify model.event.index_event_status_user_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_e_es_un';

ROLLBACK;
