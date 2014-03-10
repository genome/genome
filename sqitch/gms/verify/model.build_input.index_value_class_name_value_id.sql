-- Verify model.build_input.index_value_class_name_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_m_bi_vcn_vi';

ROLLBACK;
