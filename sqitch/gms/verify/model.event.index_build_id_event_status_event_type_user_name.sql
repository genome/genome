-- Verify model.event.index_build_id_event_status_event_type_user_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_e_bi_es_et_un';

ROLLBACK;
