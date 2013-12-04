-- Verify result.user.index_software_result_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'user_software_result_id_idx';

ROLLBACK;
