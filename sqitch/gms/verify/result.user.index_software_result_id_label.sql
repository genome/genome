-- Verify result.user.index_software_result_id_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'user_result_label_index';

ROLLBACK;
