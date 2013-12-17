-- Verify model.event.index_build_id_date_completed_event_type

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_e_bi_dc_et';

ROLLBACK;
