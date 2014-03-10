-- Verify model.build_input.index_name_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_bi_n_vi';

ROLLBACK;
