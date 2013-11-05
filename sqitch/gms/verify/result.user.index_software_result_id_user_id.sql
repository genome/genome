-- Verify result.user.index_software_result_id_user_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'user_result_id_user_id_index';

ROLLBACK;
