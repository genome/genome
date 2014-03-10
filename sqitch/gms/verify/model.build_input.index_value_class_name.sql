-- Verify model.build_input.index_value_class_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_bi_value_class_name_index';

ROLLBACK;
