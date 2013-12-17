-- Revert model.build_input.index_name_value_id

BEGIN;

DROP INDEX model.idx_m_bi_n_vi;

COMMIT;
