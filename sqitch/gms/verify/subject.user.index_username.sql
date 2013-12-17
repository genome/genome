-- Verify subject.user.index_username

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'subject_user_username_index';

ROLLBACK;
