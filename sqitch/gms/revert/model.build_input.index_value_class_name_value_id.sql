-- Revert model.build_input.index_value_class_name_value_id

BEGIN;

DROP INDEX model.idx_m_bi_vcn_vi;

COMMIT;
