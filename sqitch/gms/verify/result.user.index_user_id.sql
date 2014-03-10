-- Verify result.user.index_user_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'sru_uid_i';

ROLLBACK;
