-- Revert model.build_input.index_value_class_name

BEGIN;

DROP INDEX model.m_bi_value_class_name_index;

COMMIT;
