-- Verify model.event.index_parent_event_id_ref_seq_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_e_pei_rsi';

ROLLBACK;
